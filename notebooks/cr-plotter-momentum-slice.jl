### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ f1ee2cb0-8274-11ef-0826-f55183647219
import Pkg; Pkg.activate(Base.current_project())

# ╔═╡ 7899ae97-fbc2-43e5-ac77-c6d725f0371e
using JLD2, DataFrames

# ╔═╡ b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
using Distributions


# ╔═╡ d8a66ccd-efdc-4b85-9af3-cfa5624c88e8
# using CairoMakie
using WGLMakie

# ╔═╡ 7a050dc5-7772-4933-959f-bf4fb478fc7d
using PlutoUI

# ╔═╡ 40efcd80-db38-4db3-a193-6e65ee5c4367
using PlutoUI: Slider

# ╔═╡ 3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
using Missings

# ╔═╡ fe2b3846-c753-4685-8704-e6fb50624989
using Printf

# ╔═╡ fd47dab7-426c-44fc-8038-00e378324e41
using HypothesisTests

# ╔═╡ 49902e99-870d-4d19-afb0-1de612c185df
using StatsBase

# ╔═╡ 4ac32bef-af81-4f7e-8e97-7eac4dd2bf69
using LinearAlgebra

# ╔═╡ a5526239-2f05-4618-8868-0f552855d574
md"""
# Preamble
"""

# ╔═╡ cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
md"""
## Import packages
"""

# ╔═╡ 547aad6f-32db-405d-9886-a727f1591101
# ╠═╡ disabled = true
#=╠═╡
begin
    using AlgebraOfGraphics
    import AlgebraOfGraphics as AoG
end
  ╠═╡ =#

# ╔═╡ 5c6b130f-0a51-4131-bab7-40b059c4cc11
md"""
## Configure notebook appearance
"""

# ╔═╡ b544df91-fe2d-4396-892c-7faea2edd141
TableOfContents(depth = 6)

# ╔═╡ 4415022a-54dc-4f3d-a651-f66ae63dd051
# Increase cell width
html"""
<style>
main {
    max-width: 83%;
    padding-left: max(300px, 5%);
    padding-right: 0%;
}
</style>
"""

# ╔═╡ 8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
md"""
# Read data file
"""

# ╔═╡ d70a4da5-1589-4b41-af32-05671f27be4d
const datadir = "G:/My Drive/MC Scattering/Processed-data";

# ╔═╡ bdb9591b-b7ac-47e6-98bc-f18921bb64f9
CR_p_gdf_momentum = load_object(joinpath(datadir, "dNdp-CR-protons-momentum-split.jld2"));

# ╔═╡ 3777306e-eb41-413b-80a9-72cdc0228a94
CR_e_gdf_momentum = load_object(joinpath(datadir, "dNdp-CR-electrons-momentum-split.jld2"));

# ╔═╡ 628130bf-da25-4799-8e5e-3d2db15b1e49
md"""
# Plot Cosmic Ray data
"""

# ╔═╡ 68c8329f-501e-47df-8047-d3cbc319e705
md"""
For protons:
"""

# ╔═╡ 985a2460-3fbc-4935-af59-2e734786c973
md"""
For electrons:
"""

# ╔═╡ 59a22149-3397-4e97-9f7b-5d502aacf293
const markersize = 5;

# ╔═╡ f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
const axis_properties = (xminorgridvisible = true, yminorgridvisible = true, xlabel = "log(dN/dp)");

# ╔═╡ f86707a1-9d79-4df8-8798-3f7ea1d1797c
const bins = 90;

# ╔═╡ 377aaf8f-b909-4c42-bc77-912fd300c300
const normalization = :pdf;

# ╔═╡ 50b1a87f-49ff-4d93-aa6e-f042a87b875e
const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ 3bf64608-0fa2-4fcb-9782-fd7a8de47bda
md"""
## Sample statistics
"""

# ╔═╡ 932c2a77-0198-4df4-a4bd-30d0bda93946
md"""
### Means
"""

# ╔═╡ 5767b9ac-64c2-4d2f-ad42-961184c7edc7
md"""
### Standard deviations
"""

# ╔═╡ 7495e7e9-3d50-4401-baef-d2e3c11e6b46
md"""
### Skewness
"""

# ╔═╡ 44cb6acf-7fee-4e3e-8253-d91e5a76299a
md"""
## With Algebra of Graphics
"""

# ╔═╡ 6d5eb940-6739-4781-9dda-7433cae3cf50
Base.:*(x::Bool, l::AoG.Layer) = x ? l : AoG.zerolayer()

# ╔═╡ 2f44c2c5-7fc6-4c93-be6e-0cffb863afd4
visual_layer = AoG.histogram(Stairs; bins, normalization);

# ╔═╡ b15d71e1-ac37-422e-98f0-a0a03238fe35
map_layer = AoG.mapping(
    [:log_dNdp_cr_pf, :log_dNdp_cr_sf, :log_dNdp_cr_ISM],
    color = dims(1) => renamer(["Plasma frame", "Shock frame", "ISM frame"]));

# ╔═╡ ce8b1307-dc78-463b-9f41-04fe5dded525
md"""
## Histograms
"""

# ╔═╡ 35710ad9-f2e4-487b-be19-c29500633726
let
    idx_range = axes(CR_p_gdf_momentum,1)
    binder = @bind proton_momentum_index NumberField(idx_range, default = 13)
    min_idx, max_idx = extrema(idx_range)
    md"""
    Proton momentum slice to plot (index): $binder (min: $min_idx, max: $max_idx)
    """ # should the proton_momentum_index variable be considered a leak here?
end

# ╔═╡ d1c788a6-27ff-40d2-9bf4-1e7a4b6c48f3
df = CR_p_gdf_momentum[proton_momentum_index];

# ╔═╡ aa3a5985-3d14-4e6e-b2b2-5f7f731c3336
data = AoG.data(df);

# ╔═╡ 1333eaeb-8aae-49d5-aabc-3622b9d6ae35
layer = data * map_layer * visual_layer;

# ╔═╡ ecf80697-b786-4b02-9563-f3d082383b76
md"""
Choose which frames to plot:
- Plasma frame: $(@bind do_plot_pf CheckBox(default=true))
- Shock frame: $(@bind do_plot_sf CheckBox(default=false))
- ISM frame: $(@bind do_plot_ISM CheckBox(default=false))
"""

# ╔═╡ 7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
let
    idx_range = axes(CR_e_gdf_momentum,1)
    binder = @bind electron_momentum_index NumberField(idx_range, default = 13)
    min_idx, max_idx = extrema(idx_range)
    md"""
    Electron momentum slice to plot (index): $binder (min: $min_idx, max: $max_idx)
    """ # should the electron_momentum_index variable be considered a leak here?
end

# ╔═╡ 9ea7a3a4-987d-416d-88d1-672e3cce23c5
md"""
## dN/dp vs. iteration
"""

# ╔═╡ b7a96870-784e-4ce0-830d-d245fc16e5f4
# ╠═╡ disabled = true
#=╠═╡
let df = CR_p_gdf_momentum[proton_momentum_index]

    f = Figure()
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (protons) against iteration, momentum slice",
        xlabel = "Iteration", ylabel = "log(dN/dp)")

    do_plot_pf && scatter!(ax, df.iter, df.log_dNdp_cr_pf, label = "plasma frame"; markersize, color = color_pf_p)
    do_plot_sf && scatter!(ax, df.iter, df.log_dNdp_cr_sf, label = "shock frame"; markersize, color = color_sf_p)
    do_plot_ISM && scatter!(ax, df.iter, df.log_dNdp_cr_ISM, label = "ISM frame"; markersize, color = color_ISM_p)

    #xlims!(ax, -16, -3)
    #ylims!(ax, -100, -98)
    leg = axislegend(ax, position = :rb)
    #leg.framevisible = false
    #Legend(f[1,2], ax)
    f
end
  ╠═╡ =#

# ╔═╡ 4ac1798d-ec27-4571-9b2a-44cb432ef0d6
# ╠═╡ disabled = true
#=╠═╡
let df = CR_p_gdf_momentum[electron_momentum_index]

    f = Figure()
    ax = Axis(
        f[1,1];
        #aspect = AxisAspect(1.2),
        title = "dN/dp of Cosmic rays (electrons) against iteration, momentum slice",
        #axis_properties...,
        xlabel = "Iteration", ylabel = "log(dN/dp)")

    do_plot_pf && scatter!(ax, df.iter, df.log_dNdp_cr_pf, label = "plasma frame"; markersize, color = color_pf_e)
    do_plot_sf && scatter!(ax, df.iter, df.log_dNdp_cr_sf, label = "shock frame"; markersize, color = color_sf_e)
    do_plot_ISM && scatter!(ax, df.iter, df.log_dNdp_cr_ISM, label = "ISM frame"; markersize, color = color_ISM_e)

    #xlims!(ax, -16, -3)
    #ylims!(ax, -100, -98)
    axislegend(ax, position = :rb)
    #Legend(f[1,2], ax)
    f
end
  ╠═╡ =#

# ╔═╡ f95a0d36-5dd8-4190-98c6-06e8be2ad840
# ╠═╡ disabled = true
#=╠═╡
pcutdf = CR_p_gdf_momentum[proton_momentum_index]
  ╠═╡ =#

# ╔═╡ b88ef78f-6d6f-4b38-a9af-6da4f540f8c3
#=╠═╡
describe(pcutdf)
  ╠═╡ =#

# ╔═╡ f3132403-113d-4b30-9fd0-379d28ade3c7
md"""
# Normal distribution inference
"""

# ╔═╡ f3212b13-682f-4be4-865b-fd0f1b450aa4
Makie.update_theme!(colormap = Makie.wong_colors())

# ╔═╡ da107273-c428-4c68-80a9-8f82cb211497
md"""
# Hypothesis tests
"""

# ╔═╡ f330af91-60a6-46ac-bdc5-ec49c216fccb
md"""
Get goodness of fits
"""

# ╔═╡ 79dc57bb-d66d-4608-a775-9dfc58af1995
md"""
Plot p-values in log scale? (Uncheck for linear)
$(
    @bind plot_p_values_in_logscale CheckBox()
)
"""

# ╔═╡ 04dad413-0dc0-4ceb-81c2-e208ef082f38
p_val_yscale = plot_p_values_in_logscale ? log10 : identity;

# ╔═╡ 94a91acd-a878-4c3c-9716-8bed60bf8c6c
md"""
## Sum of squared residuals
"""

# ╔═╡ 98675d19-3b1b-4be0-9e48-ab0ffd019647
md"""
## Anderson–Darling test
"""

# ╔═╡ b499bf86-3e7a-441a-809a-934a1a8dd402
md"""
## Shapiro–Wilk test
"""

# ╔═╡ a2dca585-2b84-4958-8ba6-af51602c4d8a
sw_scores_p = let
    arr = []
    for df in CR_p_gdf_momentum
        vec = collect(skipmissing(df.log_dNdp_cr_pf))
        if length(vec) < 3
            push!(arr, missing)
            continue
        end
        score = ShapiroWilkTest(vec)
        push!(arr, score)
    end
    arr
end;

# ╔═╡ ec6883c9-b6bb-4e7e-bd6d-e65d6e06144d
sw_scores_e = let
    arr = []
    for df in CR_e_gdf_momentum
        vec = collect(skipmissing(df.log_dNdp_cr_pf))
        if length(vec) < 3
            push!(arr, missing)
            continue
        end
        score = ShapiroWilkTest(vec)
        push!(arr, score)
    end
    arr
end;

# ╔═╡ 2ab2979f-1ad4-4168-b59c-a25e57d4826a
md"""
## Kolmogorov–Smirnov test
"""

# ╔═╡ 8d03de5e-d344-4efd-b9af-dd5391028780
md"""
# Constants and functions
"""

# ╔═╡ e780481f-ffde-407f-8dff-bc289e0ceb40
"""
    fitdistribution(D, x::AbstractVector{Union{Missing,T}}) where {T}

Fit a distribution of type `D` to a vector of samples `x`.
`x` may contain `missing` values. However, if `x` contains _only_ `missing` values, `missing` is returned instead of a distribution.
This function mainly exists as a wrapper around `Distribtuions.fit` to handle `missing` values.
"""
function fitdistribution(D::Type{<:Distribution}, x::AbstractVector{Union{Missing,T}}) where {T}
    x = collect(skipmissing(x))
    if isempty(x) # don't fit to a dataset with only missings
        return missing
    end

    return Distributions.fit(D{T}, x)
end

# ╔═╡ 84d1d644-6a5b-44eb-ab4f-3b9b7171d6fe
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

# ╔═╡ e6b9701d-3d27-4c0c-b0b9-9879527f369c
normal_distrib_protons = fitdistributions(Normal, CR_p_gdf_momentum)

# ╔═╡ 2e79471f-3430-4b1c-91fe-80434de63cb2
ad_scores_p = let
    arr = []
    for (df, dist) in zip(CR_p_gdf_momentum, normal_distrib_protons.pf)
        if ismissing(dist)
            push!(arr, missing)
            continue
        end
        score = OneSampleADTest(collect(skipmissing(df.log_dNdp_cr_pf)), dist)
        push!(arr, score)
    end
    arr
end;

# ╔═╡ e8ab294c-0612-43c5-8b64-bb1ddec387ae
pf_scores = let
    arr = []
    for (df, dist) in zip(CR_p_gdf_momentum, normal_distrib_protons.pf)
        if ismissing(dist)
            push!(arr, missing)
            continue
        end
        score = ExactOneSampleKSTest(collect(skipmissing(df.log_dNdp_cr_pf)), dist)
        push!(arr, score)
    end
    arr
end

# ╔═╡ e75ea9c0-59ca-4097-b4f6-6a3af04dc308
normal_distrib_electrons = fitdistributions(Normal, CR_e_gdf_momentum)

# ╔═╡ bd8f636c-6033-434e-a220-a07397679431
ad_scores_e = let
    arr = []
    for (df, dist) in zip(CR_e_gdf_momentum, normal_distrib_electrons.pf)
        if ismissing(dist)
            push!(arr, missing)
            continue
        end
        score = OneSampleADTest(collect(skipmissing(df.log_dNdp_cr_pf)), dist)
        push!(arr, score)
    end
    arr
end;

# ╔═╡ ea3f967e-770b-4879-bddd-8d3b497344bf
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

# ╔═╡ a36ea9cf-176f-40bd-8577-cc2ea8db64af
CR_gdfstats(CR_p_gdf_momentum)

# ╔═╡ d85427f4-86ed-4c04-980a-a4152b5875e8
CR_gdfstats(CR_e_gdf_momentum)

# ╔═╡ c507291a-479c-4fa3-8956-0f841391c23f
centers(v) = (v[begin:end-1] + v[begin+1:end])/2;

# ╔═╡ 222df0cb-0760-48a2-902e-91d32e451a11
sse_scores_p = let
    arr = []
    for (df, dist) in zip(CR_p_gdf_momentum, normal_distrib_protons.pf)
        if ismissing(dist)
            push!(arr, missing)
            continue
        end
        log_dNdp = collect(skipmissing(df.log_dNdp_cr_pf))
        histogram = normalize(fit(Histogram, log_dNdp); mode=:pdf)

        x = histogram.edges |> only |> centers
        hist_y = histogram.weights
        normal_y = pdf.(dist, x)
        score = norm(hist_y - normal_y)
        push!(arr, score)
    end
    arr
end;

# ╔═╡ cbea4ff4-b132-4abb-97c6-e406a339ced6
sse_scores_e = let
    arr = []
    for (df, dist) in zip(CR_e_gdf_momentum, normal_distrib_electrons.pf)
        if ismissing(dist)
            push!(arr, missing)
            continue
        end
        log_dNdp = collect(skipmissing(df.log_dNdp_cr_pf))
        histogram = normalize(fit(Histogram, log_dNdp); mode=:pdf)

        x = histogram.edges |> only |> centers
        hist_y = histogram.weights
        normal_y = pdf.(dist, x)
        score = norm(hist_y - normal_y)
        push!(arr, score)
    end
    arr
end;

# ╔═╡ 22a91877-4339-4e6b-8f7d-0e5a80dc7c74
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

# ╔═╡ 32edc221-e586-4510-9427-977b22f62f6c
md"""
Vector of momentum slices
"""

# ╔═╡ e8406a6a-ecc2-49d2-b67a-503b4ef5764b
const proton_log_p_nat = keys(CR_p_gdf_momentum) .|> values .|> first;

# ╔═╡ 71404de8-f8b2-4d26-b7d7-41064cae1447
log_p_nat_at_slice = proton_log_p_nat[proton_momentum_index];

# ╔═╡ 4979cc00-15c1-40da-b538-021a067d1065
draw(
    layer;
    figure = (;
        title = "Histogram of protons dN/dp at p = 10^$log_p_nat_at_slice mₚc",
        titlealign = :center,
    ),
)

# ╔═╡ 89bcb29b-0b1c-4e3a-91cb-282c05df2bc5
md"""
Value of proton momentum at slice: 10^$(log_p_nat_at_slice) *m*ₚ*c*
"""

# ╔═╡ 4051e244-4c84-4983-8cb9-bc7f53daa9f6
let df = CR_p_gdf_momentum[proton_momentum_index], distribs = normal_distrib_protons
    f = Figure()
    ax = Axis(
        f[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Histogram of protons dN/dp at p = 10^$log_p_nat_at_slice mₚc",
        axis_properties...)

    if do_plot_pf
        N = df.log_dNdp_cr_pf |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "plasma frame"; bins, normalization, color = color_pf_p)

        distrib = distribs.pf[proton_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", params(distrib)...), color = color_pf_p)
        end
    end

    if do_plot_sf
        N = df.log_dNdp_cr_sf |> skipmissing |> collect
        !isempty(N) && hist!(ax, N, label = "shock frame"; bins, normalization, color = color_sf_p)
        distrib = distribs.sf[proton_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", params(distrib)...), color = color_sf_p)
        end
    end
    if do_plot_ISM
        N = df.log_dNdp_cr_ISM |> skipmissing |> collect
        !isempty(N) && hist!(ax, N, label = "ISM frame"; bins, normalization, color = color_ISM_p)
        distrib = distribs.ISM[proton_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", params(distrib)...), color = color_ISM_p)
        end
    end

    try
        axislegend(ax)
    catch e
        # axislegend has no plots to work with, because the current index doesn't have any samples. stop it complaining.
    end
    f
end

# ╔═╡ 08542eea-964a-4f1d-aae5-2b50a628588a
let
    f = Figure()
    ax = Axis(
        f[1,1];
        title = "Kolmogorov–Smirnov p-value vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)",
        yscale = p_val_yscale,
    )

    scatterlines!(ax, proton_log_p_nat, passmissing(pvalue).(pf_scores), color = color_pf_p, label = "protons, plasma frame"; markersize)
    # scatterlines!(ax, electron_log_p_nat, passmissing(pvalue).(sw_scores_e), color = color_pf_e, label = "electrons, plasma frame"; markersize)

    axislegend(ax, position = plot_p_values_in_logscale ? :cb : :lt)

    f
end

# ╔═╡ 589661b1-6a64-4db5-ac40-c1565c29c3cc
const electron_log_p_nat = keys(CR_e_gdf_momentum) .|> values .|> first;

# ╔═╡ 91bba2da-c925-4123-bb8a-c1f9be8619e9
let
    f = Figure()
    ax = Axis(
        f[1,1];
        title = "Mean vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "μ",
    )

    means = passmissing(getproperty).(normal_distrib_protons.pf, :μ)
    lines!(ax, proton_log_p_nat, means, color = color_pf_p, label = "protons, plasma frame")

    means = passmissing(getproperty).(normal_distrib_electrons.pf, :μ)
    lines!(ax, electron_log_p_nat, means, color = color_pf_e, label = "electrons, plasma frame")

    axislegend(ax)

    f
end

# ╔═╡ b6ce51e5-b4ff-49eb-83db-ecf3e8a081ac
let
    f = Figure()
    ax = Axis(
        f[1,1];
        title = "Standard deviation vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "σ",
    )
    markersize = 4

    std_devs = passmissing(getproperty).(normal_distrib_protons.pf, :σ)
    scatterlines!(ax, proton_log_p_nat, std_devs, color = color_pf_p, label = "protons, plasma frame"; markersize)

    std_devs = passmissing(getproperty).(normal_distrib_electrons.pf, :σ)
    scatterlines!(ax, electron_log_p_nat, std_devs, color = color_pf_e, label = "electrons, plasma frame"; markersize)

    axislegend(ax, position = :lt)

    f
end

# ╔═╡ adf24143-4be1-46c7-a63a-fe4dd490791d
let
    f = Figure()
    ax = Axis(
        f[1,1];
        title = "Sample skewness vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "γ",
        #yscale = log10,
    )
    # you've really gotta refactor this code
    skewness_getter(gdf) = [skewness(df[!,:log_dNdp_cr_pf] |> skipmissing |> collect) for df in gdf]

    scatterlines!(ax, proton_log_p_nat, skewness_getter(CR_p_gdf_momentum), color = color_pf_p, label = "protons, plasma frame"; markersize)
    scatterlines!(ax, electron_log_p_nat, skewness_getter(CR_e_gdf_momentum), color = color_pf_e, label = "electrons, plasma frame"; markersize)

    axislegend(ax, position = :lb)

    f
end

# ╔═╡ 6c16fc5a-7113-4b6e-abf2-de1275cceda5
log_p_nat_at_slice_e = electron_log_p_nat[electron_momentum_index];

# ╔═╡ c9b9969c-2c7f-436e-b5a1-603138a4e196
md"""
Value of electron momentum at slice: 10^$(log_p_nat_at_slice_e) *m*ₚ*c*
"""

# ╔═╡ 88822f52-aab8-4931-9091-1909da6c604b
let df = CR_e_gdf_momentum[electron_momentum_index], distribs = normal_distrib_electrons
    f = Figure()
    ax = Axis(
        f[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Histogram of electrons dN/dp at p = 10^$log_p_nat_at_slice_e mₚc",
        axis_properties...)

    if do_plot_pf
        N = df.log_dNdp_cr_pf |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "plasma frame"; bins, normalization, color = color_pf_e)
        distrib = distribs.pf[electron_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", params(distrib)...), color = color_pf_e)
        end
    end

    if do_plot_sf
        N = df.log_dNdp_cr_sf |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "shock frame"; bins, normalization, color = color_sf_e)
        distrib = distribs.sf[electron_momentum_index]
        if !ismissing(distrib)
            μ, σ = params(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", params(distrib)...), color = color_sf_e)
        end
    end
    if do_plot_ISM
        N = df.log_dNdp_cr_ISM |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "ISM frame"; bins, normalization, color = color_ISM_e)
        distrib = distribs.ISM[electron_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", params(distrib)...), color = color_ISM_e)
        end
    end

    try
        axislegend(ax, position = :lt)
    catch e
        # axislegend has no plots to work with, because the current index doesn't have any samples. stop it complaining.
    end
    f
end

# ╔═╡ e7a26d10-0e00-444d-a8f9-27874a8f821e
let
    f = Figure()
    ax = Axis(
        f[1,1];
        title = "Sum-of-Squared-Errors vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)",
        yscale = p_val_yscale,
    )

    scatterlines!(ax, proton_log_p_nat, sse_scores_p, color = color_pf_p, label = "protons, plasma frame"; markersize)
    scatterlines!(ax, electron_log_p_nat, sse_scores_e, color = color_pf_e, label = "electrons, plasma frame"; markersize)

    axislegend(ax, position = plot_p_values_in_logscale ? :cb : :ct)

    f
end

# ╔═╡ cee91c99-adc0-4185-a7c3-e2164b95a003
let
    f = Figure()
    ax = Axis(
        f[1,1];
        title = "Anderson–Darling p-value vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)",
        yscale = p_val_yscale,
    )

    scatterlines!(ax, proton_log_p_nat, passmissing(pvalue).(ad_scores_p), color = color_pf_p, label = "protons, plasma frame"; markersize)
    scatterlines!(ax, electron_log_p_nat, passmissing(pvalue).(ad_scores_e), color = color_pf_e, label = "electrons, plasma frame"; markersize)

    axislegend(ax, position = :lt)

    f
end

# ╔═╡ a49ff5ab-6077-4bb2-b694-6f3662982745
let
    f = Figure()
    ax = Axis(
        f[1,1];
        title = "Shapiro–Wilk p-value vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)",
        yscale = p_val_yscale,
    )

    scatterlines!(ax, proton_log_p_nat, passmissing(pvalue).(sw_scores_p), color = color_pf_p, label = "protons, plasma frame"; markersize)
    scatterlines!(ax, electron_log_p_nat, passmissing(pvalue).(sw_scores_e), color = color_pf_e, label = "electrons, plasma frame"; markersize)

    axislegend(ax, position = plot_p_values_in_logscale ? :cb : :lt)

    f
end

# ╔═╡ Cell order:
# ╟─a5526239-2f05-4618-8868-0f552855d574
# ╟─cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
# ╠═f1ee2cb0-8274-11ef-0826-f55183647219
# ╠═7899ae97-fbc2-43e5-ac77-c6d725f0371e
# ╠═b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
# ╠═d8a66ccd-efdc-4b85-9af3-cfa5624c88e8
# ╠═547aad6f-32db-405d-9886-a727f1591101
# ╠═7a050dc5-7772-4933-959f-bf4fb478fc7d
# ╠═40efcd80-db38-4db3-a193-6e65ee5c4367
# ╠═3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
# ╠═fe2b3846-c753-4685-8704-e6fb50624989
# ╟─5c6b130f-0a51-4131-bab7-40b059c4cc11
# ╠═b544df91-fe2d-4396-892c-7faea2edd141
# ╟─4415022a-54dc-4f3d-a651-f66ae63dd051
# ╟─8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
# ╠═d70a4da5-1589-4b41-af32-05671f27be4d
# ╠═bdb9591b-b7ac-47e6-98bc-f18921bb64f9
# ╠═3777306e-eb41-413b-80a9-72cdc0228a94
# ╟─628130bf-da25-4799-8e5e-3d2db15b1e49
# ╟─68c8329f-501e-47df-8047-d3cbc319e705
# ╟─a36ea9cf-176f-40bd-8577-cc2ea8db64af
# ╟─985a2460-3fbc-4935-af59-2e734786c973
# ╟─d85427f4-86ed-4c04-980a-a4152b5875e8
# ╠═59a22149-3397-4e97-9f7b-5d502aacf293
# ╠═f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
# ╠═f86707a1-9d79-4df8-8798-3f7ea1d1797c
# ╠═377aaf8f-b909-4c42-bc77-912fd300c300
# ╠═50b1a87f-49ff-4d93-aa6e-f042a87b875e
# ╟─3bf64608-0fa2-4fcb-9782-fd7a8de47bda
# ╟─932c2a77-0198-4df4-a4bd-30d0bda93946
# ╟─91bba2da-c925-4123-bb8a-c1f9be8619e9
# ╟─5767b9ac-64c2-4d2f-ad42-961184c7edc7
# ╟─b6ce51e5-b4ff-49eb-83db-ecf3e8a081ac
# ╟─7495e7e9-3d50-4401-baef-d2e3c11e6b46
# ╠═adf24143-4be1-46c7-a63a-fe4dd490791d
# ╟─44cb6acf-7fee-4e3e-8253-d91e5a76299a
# ╠═6d5eb940-6739-4781-9dda-7433cae3cf50
# ╠═d1c788a6-27ff-40d2-9bf4-1e7a4b6c48f3
# ╠═aa3a5985-3d14-4e6e-b2b2-5f7f731c3336
# ╠═2f44c2c5-7fc6-4c93-be6e-0cffb863afd4
# ╠═b15d71e1-ac37-422e-98f0-a0a03238fe35
# ╠═1333eaeb-8aae-49d5-aabc-3622b9d6ae35
# ╠═4979cc00-15c1-40da-b538-021a067d1065
# ╟─ce8b1307-dc78-463b-9f41-04fe5dded525
# ╠═71404de8-f8b2-4d26-b7d7-41064cae1447
# ╠═6c16fc5a-7113-4b6e-abf2-de1275cceda5
# ╟─89bcb29b-0b1c-4e3a-91cb-282c05df2bc5
# ╟─c9b9969c-2c7f-436e-b5a1-603138a4e196
# ╟─35710ad9-f2e4-487b-be19-c29500633726
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
# ╟─4051e244-4c84-4983-8cb9-bc7f53daa9f6
# ╟─7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
# ╟─88822f52-aab8-4931-9091-1909da6c604b
# ╟─9ea7a3a4-987d-416d-88d1-672e3cce23c5
# ╟─b7a96870-784e-4ce0-830d-d245fc16e5f4
# ╟─4ac1798d-ec27-4571-9b2a-44cb432ef0d6
# ╠═f95a0d36-5dd8-4190-98c6-06e8be2ad840
# ╠═b88ef78f-6d6f-4b38-a9af-6da4f540f8c3
# ╟─f3132403-113d-4b30-9fd0-379d28ade3c7
# ╠═f3212b13-682f-4be4-865b-fd0f1b450aa4
# ╠═e6b9701d-3d27-4c0c-b0b9-9879527f369c
# ╠═e75ea9c0-59ca-4097-b4f6-6a3af04dc308
# ╟─da107273-c428-4c68-80a9-8f82cb211497
# ╟─f330af91-60a6-46ac-bdc5-ec49c216fccb
# ╠═fd47dab7-426c-44fc-8038-00e378324e41
# ╟─79dc57bb-d66d-4608-a775-9dfc58af1995
# ╟─04dad413-0dc0-4ceb-81c2-e208ef082f38
# ╟─94a91acd-a878-4c3c-9716-8bed60bf8c6c
# ╠═49902e99-870d-4d19-afb0-1de612c185df
# ╠═4ac32bef-af81-4f7e-8e97-7eac4dd2bf69
# ╟─e7a26d10-0e00-444d-a8f9-27874a8f821e
# ╟─222df0cb-0760-48a2-902e-91d32e451a11
# ╟─cbea4ff4-b132-4abb-97c6-e406a339ced6
# ╟─98675d19-3b1b-4be0-9e48-ab0ffd019647
# ╟─cee91c99-adc0-4185-a7c3-e2164b95a003
# ╟─2e79471f-3430-4b1c-91fe-80434de63cb2
# ╟─bd8f636c-6033-434e-a220-a07397679431
# ╟─b499bf86-3e7a-441a-809a-934a1a8dd402
# ╟─a49ff5ab-6077-4bb2-b694-6f3662982745
# ╟─a2dca585-2b84-4958-8ba6-af51602c4d8a
# ╟─ec6883c9-b6bb-4e7e-bd6d-e65d6e06144d
# ╟─2ab2979f-1ad4-4168-b59c-a25e57d4826a
# ╟─08542eea-964a-4f1d-aae5-2b50a628588a
# ╠═e8ab294c-0612-43c5-8b64-bb1ddec387ae
# ╟─8d03de5e-d344-4efd-b9af-dd5391028780
# ╠═e780481f-ffde-407f-8dff-bc289e0ceb40
# ╠═84d1d644-6a5b-44eb-ab4f-3b9b7171d6fe
# ╠═ea3f967e-770b-4879-bddd-8d3b497344bf
# ╠═22a91877-4339-4e6b-8f7d-0e5a80dc7c74
# ╠═c507291a-479c-4fa3-8956-0f841391c23f
# ╟─32edc221-e586-4510-9427-977b22f62f6c
# ╠═e8406a6a-ecc2-49d2-b67a-503b4ef5764b
# ╠═589661b1-6a64-4db5-ac40-c1565c29c3cc
