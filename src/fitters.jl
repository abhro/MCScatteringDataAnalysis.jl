function fit_dist_to_histogram(::Type{BiNormal}, v::AbstractVector{T}; params, nbins = 150) where {T}
    v = collect(skipmissing(v))
    if isempty(v)
        return (missing, missing)
    end
    # x and y of the histogram plot if treated like a curve
    x, y = get_hist_curve(v; nbins)
    data_width = maximum(v) - minimum(v) # for setting up σ ranges

    λ_ideal, _, σ₁_ideal, μ₂_ideal, σ₂_ideal = params

    # parameter sweep
    # set up ranges for each parameter
    μ₁ = x[argmax(y)] # can get μ₁ from the peak of the histogram
    μ₂_range = sort(filter(>(38.6), x)) # consider all bin centers
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

"""
    fit(T <: Distribution, h, θ₀)

Fit a distribution given a histogram. Uses a root finder on the MLE equation.

!!! note
    This method is an instance of type piracy since `Distribution` or
    `Histogram` aren't types declared in this package.

### Arguments

- `T`: Type of distribution to fit
- `h`: Histogram
- `θ₀`: Initial guess for model parameters
"""
function StatsAPI.fit(T::Type{<:Distribution}, h::Histogram, θ₀)
    x, y = get_hist_curve(h)
    model(x, θ) = pdf(T(θ), x)

    # We want the histogram curve and the pdf to overlap as much as possible
    fit = curve_fit(model, x, y, θ₀)

    # best fit parameters
    θ_best = fit.params

    return T(θ_best)
end

"""
Given a vector `v` of sample points, fit a histogram to it, then fit a normal
distribution that histogram using least-squares optimization.
"""
function fit_dist_to_histogram(::Type{Normal}, v::AbstractVector{T}; nbins = 150) where {T}
    v = collect(skipmissing(v))
    if isempty(v)
        return missing
    end
    if length(v) == 1 # essentially a delta distribution
        return Normal(only(v), zero(T))
    end
    # x and y of the histogram plot if treated like a curve
    x, y = get_hist_curve(v; nbins)

    # Gaussian p.d.f., without creating a `Normal` because it requires σ ≥ 0.
    model(t, (μ, σ)) = @. exp(-(t - μ)^2 / 2σ^2) / (√(2π) * σ)
    # gradient of the Gaussian p.d.f. with respect to μ and σ
    function ∇model(t, (μ, σ))
        f = model(t, (μ, σ))
        J = zeros(length(t), 2)
        z = (t .- μ) / σ
        J[:, 1] .= @. f * z / σ         # ∂f/∂μ
        J[:, 2] .= @. f/σ * (z^2 - 1)   # ∂f/∂σ
        return J
    end

    # Make initial guesses as good as they can be,
    # i.e., the actual sample mean and standard deviation
    μ₀ = mean(v)
    σ₀ = std(v, mean = μ₀)

    fit = curve_fit(model, ∇model, x, y, [μ₀, σ₀])

    μ, σ = fit.param

    return Normal(μ, σ)
end

"""
    fitdistribution(DT::Type{<:Distribution}, x::AbstractVector)

Wrapper around `Distributions.fit`, but it also allows `x` to contain `missing` values.
If `x` contains _only_ missing values, or is empty, `missing` is returned.
"""
function fitdistribution(DT::Type{<:Distribution}, x::AbstractVector{Union{Missing, T}}) where {T}
    x = collect(skipmissing(x))
    if isempty(x) # don't fit to a dataset with only `missing`s
        return missing
    end

    return fit(DT{T}, x)
end

"""
    fitdistributions(fitfunc, gdf::GroupedDataFrame)

Within a `GroupedDataFrame`, call `fitdistribution` on the three columns
`:log_dNdp_cr_sf`, `:log_dNdp_cr_pf`, and `:log_dNdp_cr_ISM` and return those distributions.

### Returns
- A 3-element `NamedTuple` containing:
  - `sf`: The distributions found by fitting to the `log_dNdp_cr_sf` column in each group.
  - `pf`: The distributions found by fitting to the `log_dNdp_cr_pf` column in each group.
  - `ISM`: The distributions found by fitting to the `log_dNdp_cr_ISM` column in each group.
"""
function fitdistributions(fitfunc, gdf::GroupedDataFrame)

    sf = Vector{Any}(undef, length(gdf))
    pf = Vector{Any}(undef, length(gdf))
    ISM = Vector{Any}(undef, length(gdf))

    for (i, df) in enumerate(gdf)
        # fit a distribution to the {shock,plasma,ISM} frame data
        cursf = fitfunc(df.log_dNdp_cr_sf)
        curpf = fitfunc(df.log_dNdp_cr_pf)
        curISM = fitfunc(df.log_dNdp_cr_ISM)

        sf[i] = cursf
        pf[i] = curpf
        ISM[i] = curISM
    end

    # narrow the element type of each vector of distributions
    sf = Vector{Union{Set(typeof.(sf))...}}(sf)
    pf = Vector{Union{Set(typeof.(pf))...}}(pf)
    ISM = Vector{Union{Set(typeof.(ISM))...}}(ISM)

    return (; sf, pf, ISM)
end

"""
    fitnormal(x::AbstractVector; corrected = false)

Analogue of `fitdistribution`, but directly constructs a `Normal` using the mean and variance.

Note: Does not use the Bessel correction for the variance by default (set `corrected`)
"""
function fitnormal(x::AbstractVector; corrected = false)
    x = collect(skipmissing(x))
    if isempty(x) # don't fit to a dataset with only missing values
        return missing
    end

    μ = mean(x)
    σ = std(x; corrected, mean = μ)
    if isnan(μ) || isnan(σ)
        return missing
    end
    return Normal(μ, σ)
end

"""
    specific_width_histogram_fits(gdf, width, col = :log_dNdp_cr_pf, normalization = :pdf)

Like `fit(Histogram, v, ...)`, but specify bin `width` instead of `nbins`.
"""
function specific_width_histogram_fits(gdf, width, col = :log_dNdp_cr_pf, normalization = :pdf)
    hists = Vector{Any}(undef, length(gdf))
    for (i, df) in enumerate(gdf)
        v = df[!, col] |> skipmissing |> collect
        if length(v) < 3
            hists[i] = missing
            continue
        end
        hist = fit(Histogram, v, edges(v, width))
        hist = normalize(hist; mode = normalization)
        hists[i] = hist
    end
    return hists
end

"""
    edges(v, width)

Given a sample vector `v` and a bin width `width`, return a range spanning all
the values of `v` where adjacent elements are `width` apart.

TODO: write better docstring
"""
function edges(v, width)
    xmin, xmax = extrema(skipmissing(v))
    return range(xmin, xmax, step = width)
end
