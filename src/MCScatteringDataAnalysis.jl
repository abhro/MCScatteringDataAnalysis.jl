module MCScatteringDataAnalysis

using LinearAlgebra
using StatsBase
using Distributions
using DataFrames
using BiNormalDistributions

sse(ŷ, y) = sum((y - ŷ).^2)
export sse

function fit_dist_to_histogram(v::AbstractVector{T}; params, nbins = 150) where T
    h = normalize(StatsBase.fit(StatsBase.Histogram, v; nbins), mode = :pdf)
    edges = only(h.edges)
    # x and y of the histogram plot if treated like a curve
    x = centers(edges)
    #x = edges[begin+1:end]
    y = h.weights
    data_width = maximum(v) - minimum(v) # for setting up sigma ranges
    @debug("Calculated histogram", h, edges, x, y, data_width)

    lambda_ideal, _, sigma_1_ideal, mu_2_ideal, sigma_2_ideal = params

    # parameter sweep
    # set up ranges for each parameter
    μ₁ = x[argmax(y)] # can get μ₁ from the peak of the histogram
    μ₂_range = sort(filter(>(38.6), vcat(edges, x))) # consider all bin edges and centers
    #σ₁_range = range(0, data_width/2, length = 50)
    σ₁_range = range(0.1, 0.2, step = 0.01)
    #σ₂_range = σ₁_range # use same range for both
    σ₂_range = range(0.7, 0.8, step = 0.01)
    λ_range = 0.5:0.01:1.0
    @debug("Using ranges", μ₁, μ₂_range, σ₁_range, σ₂_range, λ_range)

    # starting guess: single Gaussian, unit variance
    local best_model = BiNormal(1.0, μ₁, one(T), zero(T), one(T))
    local best_fit_score = Inf

    # for each parameter (Cartesian product of all ranges)
    for μ₂ in μ₂_range, σ₁ in σ₁_range, σ₂ in σ₂_range, λ in λ_range
        # create the model
        model = BiNormal(λ, μ₁, σ₁, μ₂, σ₂)
        # calculate the ŷ produced by the new parameter set
        ŷ = pdf.(model, x)
        # calculate the goodness of fit
        fit_score = sse(ŷ, y)

        # if better than current best fit, save model and new fitness score
        if fit_score < best_fit_score
            best_model = model
            best_fit_score = fit_score

            #@debug("Found new best model", model, fit_score)
        end
    end
    return (best_model, best_fit_score)
end
export fit_dist_to_histogram

"""
    fitdistribution(DT::Type{<:Distribution}, x::AbstractVector)

Wrapper around `Distributions.fit`, but it also allows `x` to contain `missing` values.
If `x` contains _only_ missing values, or is empty, `missing` is returned.
"""
function fitdistribution(DT::Type{<:Distribution}, x::AbstractVector{Union{Missing,T}}) where {T}
    x = collect(skipmissing(x))
    if isempty(x) # don't fit to a dataset with only `missing`s
        return missing
    end

    return Distributions.fit(DT{T}, x)
end
export fitdistribution

"""
    fitdistributions(DT::Type{<:Distribution}, gdf::GroupedDataFrame)

Within a `GroupedDataFrame`, call `fitdistribution` on the three columns
`:log_dNdp_cr_sf`, `:log_dNdp_cr_pf`, and `:log_dNdp_cr_ISM` and return those distributions.

### Returns
- A 3-element `NamedTuple` containing:
  - `sf`: The distributions found by fitting to the `log_dNdp_cr_sf` column in each group.
  - `pf`: The distributions found by fitting to the `log_dNdp_cr_pf` column in each group.
  - `ISM`: The distributions found by fitting to the `log_dNdp_cr_ISM` column in each group.
"""
function fitdistributions(DT::Type{<:Distribution}, gdf::GroupedDataFrame)

    DistArrayType = Vector{Union{Missing,Nothing,Distribution}}

    sf = DistArrayType(undef, length(gdf))
    pf = DistArrayType(undef, length(gdf))
    ISM = DistArrayType(undef, length(gdf))

    for (i, df) in enumerate(gdf)
        # fit a distribution to the shock frame data
        cursf = fitdistribution(DT, df.log_dNdp_cr_sf)
        # fit a distribution to the plasma frame data
        curpf = fitdistribution(DT, df.log_dNdp_cr_pf)
        # fit a distribution to the ISM frame data
        curISM = fitdistribution(DT, df.log_dNdp_cr_ISM)

        sf[i] = cursf
        pf[i] = curpf
        ISM[i] = curISM
    end

    (; sf, pf, ISM)
end
export fitdistributions

"""
    CR_gdfstats(gdf)

For a `GroupedDataFrame` of dN/dp values, compute various statistics
for each frame, grouped by momentum.
"""
function CR_gdfstats(gdf)
    n = length(gdf)
    log_p = zeros(n)
    nrows = zeros(Int, n)
    n_pf_samples = zeros(Int, n)
    n_sf_samples = zeros(Int, n)
    n_ISM_samples = zeros(Int, n)

    for (i, df) in enumerate(gdf)
        log_p[i] = keys(gdf)[i] |> values |> first
        nrows[i] = nrow(df)
        n_pf_samples[i] = count(!ismissing, df.log_dNdp_cr_pf)
        n_sf_samples[i] = count(!ismissing, df.log_dNdp_cr_sf)
        n_ISM_samples[i] = count(!ismissing, df.log_dNdp_cr_ISM)
    end
    return DataFrame(;
        log_p,
        nrows,
        n_pf_samples,
        n_sf_samples,
        n_ISM_samples,
    )
end
export CR_gdfstats

"""
    SSE_hist(occurrences, dist; relative = true, bias)

Given a list of `occurrences`, and a pre-fit distribution `dist`,
compute the L2-norm of errors of the histogram and the pdf of `dist`.
The `pdf` is evaluated at the center of each bin.

### Arguments
- `occurrences`: Array of data points
- `dist`: Distribution to compare to. The distribution is treated as the ground
  truth since the histogram from `occurrences` can have bins with 0 density.
- `relative`: Whether the relative (or absolute) error should be calculated
- `bias`: Offset term used when calculating the relative error to prevent division by
  small floating point numbers. Defaults to `eps(T)` where `T` is the type of a data point.
"""
function SSE_hist(occurrences, dist; relative = true, bias = eps(eltype(occurrences)))
    occurrences = collect(skipmissing(occurrences))
    histogram = normalize(fit(Histogram, occurrences); mode=:pdf)

    x = histogram.edges |> only |> centers
    hist_y = histogram.weights
    dist_y = pdf.(dist, x)
    residuals = hist_y - dist_y
    if relative
        residuals ./= dist_y .+ bias
    end
    score = norm(residuals)
    return score
end
export SSE_hist

"""
    centers(v)

Return array of each element being the center (mean) of adjacent elements
"""
centers(v) = (v[begin:end-1] + v[begin+1:end])/2;
export centers

"""
    get_hist_curve(occurrences; nbins)

Get a pdf normalized histogram with `nbins` bins, treated as a curve with the
x-values as the bin centers, and the y-values as the value of the pdf at the bin center.

### Returns
- `x`: center of bins
- `y`: pdf at each `x`.
"""
function get_hist_curve(occurrences; nbins)
    histogram = normalize(fit(Histogram, occurrences; nbins); mode=:pdf)

    x = histogram.edges |> only |> centers
    hist_y = histogram.weights
    return x, hist_y
end
export get_hist_curve

"""
    fitnormal(x::AbstractVector)

Analogue of `fitdistribution`, but directly constructs a `Normal` using the mean and variance.
"""
function fitnormal(x::AbstractVector{Union{Missing,T}}) where {T}
    x = collect(skipmissing(x))
    if isempty(x) # don't fit to a dataset with only missings
        return missing
    end

    μ = mean(x)
    σ = std(x, corrected = false, mean = μ)
    if isnan(μ) || isnan(σ)
        return missing
    end
    return Normal(μ, σ)
end
export fitnormal

"""
    fitnormals(gdf::GroupedDataFrame)

Analogue of `fitdistributions`, but directly constructs a `Normal` using the mean and variance.
"""
function fitnormals(gdf::GroupedDataFrame)

    DistArrayType = Vector{Union{Missing,Nothing,Normal}}

    sf = DistArrayType(undef, length(gdf))
    pf = DistArrayType(undef, length(gdf))
    ISM = DistArrayType(undef, length(gdf))

    for (i, df) in enumerate(gdf)
        # fit a distribution to the shock frame data
        cursf = fitnormal(df.log_dNdp_cr_sf)
        # fit a distribution to the plasma frame data
        curpf = fitnormal(df.log_dNdp_cr_pf)
        # fit a distribution to the ISM frame data
        curISM = fitnormal(df.log_dNdp_cr_ISM)

        sf[i] = cursf
        pf[i] = curpf
        ISM[i] = curISM
    end

    (; sf, pf, ISM)
end
export fitnormals
end
