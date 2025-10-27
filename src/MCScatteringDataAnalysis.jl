module MCScatteringDataAnalysis

using LinearAlgebra
using StatsBase
using Distributions
using DataFrames: DataFrame, GroupedDataFrame, nrow
using BiNormalDistributions: BiNormal
using HypothesisTests: OneSampleADTest, ExactOneSampleKSTest, ShapiroWilkTest
using LsqFit: curve_fit

"""
    ColumnSpecification(colname, eltype, uses_sentinels = false)

Helper datatype for parsing text data created by the mc\\_cr program.
"""
@kwdef struct ColumnSpecification
    colname::Symbol
    eltype::Type
    uses_sentinels::Bool
    ColumnSpecification(colname, eltype, uses_sentinels = false) =
        new(colname, eltype, uses_sentinels)
end
const CS = ColumnSpecification
colname(cs::CS) = cs.colname
name(cs::CS) = cs.colname
Base.eltype(cs::CS) = cs.eltype
uses_sentinels(cs::CS) = cs.uses_sentinels

include("colspecs.jl")

sse(ŷ, y) = sum((y - ŷ).^2)
export sse

include("fitters.jl")
export fit_dist_to_histogram
export fitdistribution
export fitdistributions
export fitnormal

include("bhattacharya_distance.jl")
export bcdistance

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
function SSE_hist(occurrences::AbstractVector{T}, dist;
                  nbins = 90, relative = true, bias = eps(T)) where {T}
    occurrences = collect(skipmissing(occurrences))
    x, hist_y = get_hist_curve(occurrences; nbins)
    dist_y = pdf.(dist, x)
    residuals = hist_y - dist_y
    if relative
        #residuals ./= dist_y .+ bias
        residuals ./= max.(dist_y, hist_y)
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
    histogram = normalize(StatsBase.fit(Histogram, occurrences; nbins); mode=:pdf)

    x = histogram.edges |> only |> centers
    hist_y = histogram.weights
    return x, hist_y
end
export get_hist_curve

"""
    get_sse_scores(gdf, dists)

Get Root-sum-square-error scores for a `GroupedDataFrame` `gdf` when compared
with a list of distributions `dists`.
"""
get_sse_scores(gdf, dists; col) = get_onesample_scores(SSE_hist, gdf, dists; col)
export get_sse_scores

"""
    get_sw_scores(gdf)

Get Shapiro―Wilk scores for a `GroupedDataFrame` `gdf`.
"""
function get_sw_scores(gdf; col)
    n = length(gdf)
    arr = Vector{Union{ShapiroWilkTest,Missing}}(undef, n)
    for (i, df) in enumerate(gdf)
        vec = df[!, col] |> skipmissing |> collect
        if length(vec) < 3
            arr[i] = missing
            continue
        end
        score = ShapiroWilkTest(vec)
        arr[i] = score
    end
    return arr
end
export get_sw_scores

"""
    get_ad_scores(gdf, dists)

Get Anderson―Darling scores for a `GroupedDataFrame` `gdf` when compared with a
list of distributions `dists`.
"""
get_ad_scores(gdf, dists; col) = get_onesample_scores(OneSampleADTest, gdf, dists; col)
export get_ad_scores

"""
    get_ks_scores(gdf, dists; col)

Get Kolmogorov―Smirnov scores for a `GroupedDataFrame` `gdf` when compared with a
list of distributions `dists`.
"""
get_ks_scores(gdf, dists; col) = get_onesample_scores(ExactOneSampleKSTest, gdf, dists; col)
export get_ks_scores

function get_onesample_scores(test, gdf, dists; col)
    n = length(gdf)
    arr = Vector{Any}(undef, n)
    for (i, (df, dist)) in enumerate(zip(gdf, dists))
        if ismissing(dist)
            arr[i] = missing
            continue
        end
        vec = df[!, col] |> skipmissing |> collect
        score = test(vec, dist)
        arr[i] = score
    end
    arr = Vector{Union{Set(typeof.(arr))...}}(arr)
    return arr
end
end
