module MCScatteringDataAnalysis

export CR_gdfstats, SSE_hist, fit_dist_to_histogram, fitdistribution, fitdistributions


function fit_dist_to_histogram(v::AbstractVector{T}; nbins = 150) where T
    h = normalize(StatsBase.fit(StatsBase.Histogram, v; nbins), mode = :pdf)
    edges = only(h.edges)
    # x and y of the histogram plot if treated like a curve
    x = centers(edges)
    #x = edges[begin+1:end]
    y = h.weights
    data_width = maximum(v) - minimum(v) # for setting up sigma ranges
    @debug("Calculated histogram", h, edges, x, y, data_width)

    lambda_ideal, _, sigma_1_ideal, mu_2_ideal, sigma_2_ideal =params(manual_bn)

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

    # starting guess: single gaussian, unit variance
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

function fitdistribution(DT::Type{<:Distribution}, x::AbstractVector{Union{Missing,T}}) where {T}
    x = collect(skipmissing(x))
    if isempty(x) # don't fit to a dataset with only `missing`s
        return missing
    end

    return Distributions.fit(DT{T}, x)
end

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

"""
    CR_gdfstats(gdf)

For a `GroupedDataFrame` of dN/dp values, compute various statistics grouped by momentum.
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

"""
    SSE_hist(occurrences, dist)

Given a list of `occurences`, and a pre-fit distribution `dist`,
compute the L2-norm of errors of the histogram and the pdf of `dist`.
The `pdf` is evaluated at the center of each bin.
"""
function SSE_hist(occurrences, dist)
    occurrences = collect(skipmissing(occurrences))
    histogram = normalize(fit(Histogram, occurrences); mode=:pdf)

    x = histogram.edges |> only |> centers
    hist_y = histogram.weights
    dist_y = pdf.(dist, x)
    score = norm(hist_y - dist_y)
    return score
end
end
