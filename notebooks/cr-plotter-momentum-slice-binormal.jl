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

# ╔═╡ 325d0adb-9001-4f1c-98fb-cfce042a09ca
using MCScatteringDataAnalysis

# ╔═╡ 1e0808e1-a106-4f0e-8649-13989b8ca855
using CairoMakie
#using WGLMakie

# ╔═╡ 547aad6f-32db-405d-9886-a727f1591101
begin
    using AlgebraOfGraphics
    import AlgebraOfGraphics as AoG
end

# ╔═╡ 7a050dc5-7772-4933-959f-bf4fb478fc7d
using PlutoUI

# ╔═╡ 40efcd80-db38-4db3-a193-6e65ee5c4367
using PlutoUI: Slider

# ╔═╡ 3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
using Missings

# ╔═╡ fe2b3846-c753-4685-8704-e6fb50624989
using Printf

# ╔═╡ 75e49b40-bff0-48f5-ab57-28b185f63cc9
using BiNormalDistributions

# ╔═╡ de7bbc47-64c6-4e3c-bd2e-71f189225d52
using KernelDensity

# ╔═╡ f297f99f-5f0a-4b7b-9a8b-3d7cbf37e102
using Peaks

# ╔═╡ 5445b481-0ea6-4543-b554-7ccd257cbac6
using NonlinearSolve, Optimization

# ╔═╡ 66a40e6b-a7f6-4b78-b5e2-f343fb02f4fe
using StaticArrays

# ╔═╡ 77733d10-ba25-4de5-b8a3-c52a4130227c
using LinearAlgebra

# ╔═╡ 3b6808a5-d4fc-4347-91ec-3cd389b13534
using Optim

# ╔═╡ fae42fad-f89b-4527-8bc2-8747951ea405
using NLsolve

# ╔═╡ ff76d92a-c68b-450e-b742-9b1706c5d310
using OptimizationOptimisers

# ╔═╡ f12b556e-5b3b-47f8-a7e6-4547e6c13d39
using ExpectationMaximization

# ╔═╡ 4153a601-06c3-4126-ace6-d354064e03f5
using Random

# ╔═╡ 4a0e2184-0950-4b19-9b8b-061150d17ec5
md"""
# Plot fluxes for each momentum slice with binormal distribution estimates
"""

# ╔═╡ a5526239-2f05-4618-8868-0f552855d574
md"""
## Preamble
"""

# ╔═╡ cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
md"""
### Import packages
"""

# ╔═╡ 334b4ffc-1c5d-4743-88fb-ab383a3e6f80
md"""
### Configure notebook appearance
"""

# ╔═╡ b544df91-fe2d-4396-892c-7faea2edd141
TableOfContents(depth = 6)

# ╔═╡ 4415022a-54dc-4f3d-a651-f66ae63dd051
# Increase cell width
html"""<style>
main {
    max-width: 83%;
    padding-left: max(300px, 5%);
    padding-right: 0%;
}
</style>"""

# ╔═╡ 8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
md"""
## Read data file
"""

# ╔═╡ c159f801-b129-4919-85ef-29eedf977f14
const datadir = "G:/My Drive/MC Scattering/Processed-data";

# ╔═╡ bdb9591b-b7ac-47e6-98bc-f18921bb64f9
CR_p_gdf_momentum = load_object(joinpath(datadir, "dNdp-CR-protons-momentum-split.jld2"));

# ╔═╡ 3777306e-eb41-413b-80a9-72cdc0228a94
CR_e_gdf_momentum = load_object(joinpath(datadir, "dNdp-CR-electrons-momentum-split.jld2"));

# ╔═╡ 628130bf-da25-4799-8e5e-3d2db15b1e49
md"""
## Plot Cosmic Ray data
"""

# ╔═╡ bfc6a515-8189-487b-be08-746d865a78ae
md"""
For protons
"""

# ╔═╡ c5947192-0fa5-4063-8af2-74febf514b8b
CR_gdfstats(CR_p_gdf_momentum)

# ╔═╡ ea647872-9dc3-4fb9-9499-e396127703b2
md"""
For electrons:
"""

# ╔═╡ 60ee4f38-e85f-4a2d-b17c-579531588058
CR_gdfstats(CR_e_gdf_momentum)

# ╔═╡ ecf80697-b786-4b02-9563-f3d082383b76
md"""
Choose which frames to plot:
- Plasma frame: $(@bind do_plot_pf CheckBox(default=true))
- Shock frame: $(@bind do_plot_sf CheckBox(default=false))
- ISM frame: $(@bind do_plot_ISM CheckBox(default=false))
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
## Bi-normal distribution inference
"""

# ╔═╡ 90e850f3-7b48-441a-92ab-1e1f6bf04e9a
md"""
Pick a specific momentum index to work with
"""

# ╔═╡ 18ae83a7-98e2-4ef0-b21c-cac428146188
testset_index = 67;

# ╔═╡ e0cc631f-28a1-42db-84fc-9e7dcc9387bf
CR_p_gdf_momentum[testset_index]

# ╔═╡ b63ff630-624b-4e9b-bc03-dc32fd691b05
testset = CR_p_gdf_momentum[testset_index].log_dNdp_cr_pf |> skipmissing|> collect

# ╔═╡ 96b36184-f98e-4b31-a2d5-1754bb40d84a
filter(:log_dNdp_cr_pf => >(38.7), CR_p_gdf_momentum[testset_index])

# ╔═╡ 3596bac9-5797-40c6-a4da-cdcc1cc9a451
testset

# ╔═╡ 7e776afd-fb12-4c90-a489-966541540599
n2_tentative = Distributions.fit(Normal, filter(>(38.7), testset))

# ╔═╡ ac9a9859-cf1f-4084-b784-47315c0e18c1
bn_tentative = BiNormal(
    0.99,
    Distributions.fit(Normal, testset),
    n2_tentative)

# ╔═╡ e3dab56a-5560-4ea8-84d7-53ce88cedc1c
md"""
### Manual fit
"""

# ╔═╡ d21499d1-f010-444a-96c2-1dee378496e7
manual_bn = BiNormal(0.995, 38.212839223905315, 0.17, 38.95, 0.05)

# ╔═╡ 247d55e1-c7a4-4ccc-bce3-694f4e46dc14
md"""
SSE of manual fit
"""

# ╔═╡ 7d2a121d-7e87-471b-8cc0-034796151b84
md"""
### Parameter sweep on histogram discretization
"""

# ╔═╡ 9bc9c0f1-c3be-4944-9fec-575c4fda3ce5
import StatsBase

# ╔═╡ 4eeb0f0b-311d-410b-b528-cbcb6f7490a7
q = normalize(StatsBase.fit(StatsBase.Histogram, testset, nbins = 90), mode = :pdf)

# ╔═╡ 860ad43f-8683-481e-b0d2-06194ebc1af9
let
    ŷ = pdf.(manual_bn, centers(q.edges |> first))
    y = q.weights
    sse(ŷ, y)
end

# ╔═╡ bc44add8-d20f-4e67-ae68-7af945020d55
argmax(q.weights)

# ╔═╡ 6eeb6453-dee0-45eb-89ec-19b7cb2d26c1
plot(q)

# ╔═╡ bd7fa9ce-3049-44d3-844f-df048003bfc5
q.weights |> Print

# ╔═╡ 53fbfeb5-7993-453f-bd57-2d3c409ed46e
q.edges

# ╔═╡ e6a6cbe6-d8e4-40e1-8c29-aea7703f35a9
brute_fitted = fit_dist_to_histogram(testset, params = params(manual_bn))

# ╔═╡ a97412e6-9681-4afa-8ceb-6f37f2f6dd0b
brute_fit_dist = first(brute_fitted)

# ╔═╡ d72cc184-01ce-440d-90e3-6977f9b8af7e
centers(q.edges |> only)

# ╔═╡ df94b5b9-959c-49b8-b0b0-d8c965c61a9b
q.weights

# ╔═╡ a10b1208-3349-4a66-abc6-097bd0d5acd4
import LsqFit

# ╔═╡ 82855970-94bf-4a67-8504-f9536f03722c
md"""
domain transformation:
```math
β = \ln \frac{2λ-1}{2-2λ},
\quad
λ = \frac{2+e^{-β}}{2+2e^{-β}}
```

Takes constraints from ``λ ∈ [1/2, 1]`` to ``β ∈ ℝ``.
"""

# ╔═╡ c1b47c6c-66ea-4014-9c04-8aa142449178
md"""
### Mode finding through kernel density estimates
"""

# ╔═╡ 8999b23b-4357-4655-baa6-273b218006b7
CairoMakie.density(testset)

# ╔═╡ 43a719e7-97d5-4e36-ba6f-3a4ec4b02463
##density_maxes = findmaxima(testset_kde.density) |> peakproms

# ╔═╡ b822d103-e6c1-4b76-86ea-84eb84736133
testset_kde, density_maxes = BiNormalDistributions.kdemaxes(testset, 2)

# ╔═╡ dbcc47ad-952e-475b-9657-f7fd280de743
let
    f, ax, _ = plot(testset_kde, label = "KDE")
    plot!(ax, testset_kde.x[density_maxes.indices], density_maxes.heights, label = "modes")
    axislegend(ax)
    f
end

# ╔═╡ d6516ed8-0a21-4509-a1b4-34f6521ab222
testset_kde.x[density_maxes.indices]

# ╔═╡ 5252c2c6-969d-45c1-839c-32db557aa4b8
md"""
### Fit Kernel density estimate curve
"""

# ╔═╡ 5b7baac9-e657-4666-be4d-62233362aa09
testset_kde.x |> length

# ╔═╡ 54c75433-c9e3-4d13-863e-3a1aa51f5e3e
testset_kde.density |> length

# ╔═╡ 2f5f0461-0101-41fb-b785-4cd96d455476


# ╔═╡ 04913ba6-6fa8-4c68-aeb3-5d9f189894ee
md"""
Some imports for some reason
"""

# ╔═╡ b9892967-520c-40e6-8cf1-3b1eb081ce04
md"""
### Log-likelihood maximization
"""

# ╔═╡ a08fe436-01cf-498e-8976-6e2c3173ca11
md"""
### Expectation maximization
"""

# ╔═╡ 1360285b-8b6a-4d1c-bbb5-c6acfeddb8b6
testset_mle_fit_distrib = fit_mle(MixtureModel([Normal(), Normal()], ), testset; method = StochasticEM(MersenneTwister(2)))

# ╔═╡ 788836a8-e168-4eed-b5cd-3522e43b80a6
let f = Figure()
    ax = Axis(f[1,1], xminorgridvisible = true, yminorgridvisible = true)
    ##stephist!(ax, testset; bins, normalization = :pdf, label = "stephist")
    plot!(ax, testset_kde, label = "KDE")
    xplt = range(extrema(testset)..., length = 1000)
    lines!(ax, xplt, testset_mle_fit_distrib, label = "EM fit")
    ##lines!(ax, xplt, pdf.(mixture_model_test, xplt))
    axislegend(ax)
    f
end

# ╔═╡ 8d03de5e-d344-4efd-b9af-dd5391028780
md"""
## Constants and functions
"""

# ╔═╡ 6d5eb940-6739-4781-9dda-7433cae3cf50
Base.:*(x::Bool, l::AoG.Layer) = x ? l : AoG.zerolayer()

# ╔═╡ 0aab1add-5285-4da7-b4eb-d1445b96b035
brute_fit_dist.N₁ * brute_fit_dist.λ

# ╔═╡ 7f0c2a67-1631-47d7-9d81-f87b44eab1c4
function modelfunc(x, (β, μ₁, s₁, μ₂, s₂))
    λ = (2 + exp(-β)) / (2 + 2exp(-β))
    return pdf(BiNormal(λ, μ₁, s₁^2, μ₂, s₂^2), x)
end

# ╔═╡ 00f44bb7-7074-457c-aef8-566da755d748
hist_curve_fit = let
    x_data = centers(q.edges |> only)
    LsqFit.curve_fit(
        modelfunc, x_data, q.weights,
        [1.0, mean(x_data), 1.0, mean(x_data), 1.0])
end

# ╔═╡ 7140dc51-55ca-437a-a6da-8812ffc35332
kde_curve_fit = let
    x_data = testset_kde.x
    LsqFit.curve_fit(
        modelfunc, testset_kde.x, testset_kde.density,
        [1.0, mean(x_data), 1.0, mean(x_data), 1.0])
end

# ╔═╡ c8a9815e-ed1b-44de-8d8b-9aee518cfe4e
hist_curve_fit_distrib = let
    β, μ₁, s₁, μ₂, s₂ = hist_curve_fit.param
    λ = (2 + exp(-β)) / (2 + 2exp(-β))
    BiNormal(λ, μ₁, s₁^2, μ₂, s₂^2)
end

# ╔═╡ 78e6146f-74f1-4033-9a18-1e5f8cfdd8cd
kde_curve_fit_distrib = let
    β, μ₁, s₁, μ₂, s₂ = kde_curve_fit.param
    λ = (2 + exp(-β)) / (2 + 2exp(-β))
    BiNormal(λ, μ₁, s₁^2, μ₂, s₂^2)
end

# ╔═╡ 8572d3a4-405c-438c-9dfc-0d37222eee9b
let f = Figure()
    ax = Axis(f[1,1], xminorgridvisible = true, yminorgridvisible = true)
    ##stephist!(ax, testset; bins, normalization = :pdf)
    plot!(ax, testset_kde, label = "kde")
    xplt = range(extrema(testset)..., length = 1000)
    lines!(ax, xplt, kde_curve_fit_distrib, label = "Curve fit on kde, LsqFit.jl")
    #lines!(ax, xplt, pdf.(mixture_model_test, xplt))
    axislegend(ax)
    f
end

# ╔═╡ 2f36a2d9-d9f5-4bf6-8fa4-c1a07532e8bb
const proton_indices = axes(CR_p_gdf_momentum, 1);

# ╔═╡ baf25a3a-0d13-409e-b5d9-5a1171da28b2
const electron_indices = axes(CR_e_gdf_momentum, 1);

# ╔═╡ 4e465b9b-b2a1-42c0-ab78-ea4f620dbe30
const proton_index_binder = @bind proton_momentum_index NumberField(proton_indices, default = 13);

# ╔═╡ 71404de8-f8b2-4d26-b7d7-41064cae1447
log_p_nat_at_slice_p = keys(CR_p_gdf_momentum)[proton_momentum_index] |> values |> only;

# ╔═╡ 35710ad9-f2e4-487b-be19-c29500633726
md"""
Proton momentum slice to plot (index): $(proton_index_binder) (min: $(minimum(proton_indices)), max: $(maximum(proton_indices)))

Value of proton momentum at slice: 10^$(log_p_nat_at_slice_p) _m_ₚ_c_
"""

# ╔═╡ c21810ac-c7d7-4faf-8b2d-8985adb268da
const electron_index_binder = @bind electron_momentum_index NumberField(electron_indices, default = 13);

# ╔═╡ cef8f0a4-0967-4e86-bfde-7fa84c474e31
log_p_nat_at_slice_e = keys(CR_p_gdf_momentum)[electron_momentum_index] |> values |> only;

# ╔═╡ 7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
md"""
Electron momentum slice to plot (index): $electron_index_binder (min: $(minimum(electron_indices)), max: $(maximum(electron_indices)))

Value of electron momentum at slice: 10^$(log_p_nat_at_slice_e) *m*ₚ*c*
"""

# ╔═╡ 59a22149-3397-4e97-9f7b-5d502aacf293
const markersize = 6;

# ╔═╡ f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
const axis_properties = (xminorgridvisible = true, yminorgridvisible = true, xlabel = "log(dN/dp)");

# ╔═╡ 4d245ac7-6329-457a-970e-8a8aa23775dc
#let
with_theme() do
    f = Figure()
    ax = Axis(f[1,1]; axis_properties...)
    x = centers(only(q.edges))
    lines!(ax, x, q.weights, label = "Histogram*, with bin centers")
    modevalue, modeidx = findmax(q.weights)
    @info "Got values" modevalue modeidx x[modeidx]
    scatter!(ax, x[modeidx], modevalue, label = "Peak")
    axislegend(ax)
    f
end

# ╔═╡ f86707a1-9d79-4df8-8798-3f7ea1d1797c
const bins = 90;

# ╔═╡ cf870504-0f29-4354-9a4a-76971459aeba
let testset = filter(>(38.7), testset)
    f = Figure()
    ax = Axis(f[1,1]; axis_properties...)
    hist!(ax, testset; bins, normalization = :pdf)
    #hist(; bins, normalization = :pdf)
    #xplt = range(extrema(testset)..., length = 1000)
    #lines!(ax, xplt, pdf.(mixture_model_test, xplt) * 60)
    f
end

# ╔═╡ 0b1c1d4f-6ffd-423b-bf9b-31b229488038
with_theme() do
    f = Figure()
    ax = Axis(f[1,1], title = "Fit after filtering out main distrib and manually adjusting λ"; axis_properties...)
    stephist!(ax, testset, normalization = :pdf; bins, label = "Test set")
    xs = range(extrema(testset)..., length=1000)
    #plot!(ax, xs, bn_tentative, label = "semi-Manual BiNormal")
    plot!(ax, xs, manual_bn, label = "Manual BiNormal", linewidth=1, color = :orange)
    axislegend(ax)
    f
end

# ╔═╡ c305f828-96c5-4839-9524-6a890a5d68fa
#let
with_theme(Makie.theme(nothing)) do
    f = Figure()
    ax = Axis(f[1,1]; axis_properties...)
    stephist!(ax, testset; bins, normalization = :pdf)
    xplt = range(extrema(testset)..., length = 1000)
    λ = brute_fit_dist.λ
    @info λ
    lines!(ax, xplt, brute_fit_dist, label = "Curve fit on histogram, parameter sweep")
    ##lines!(ax, xplt, λ * pdf.(brute_fit_dist.N₁, xplt), label = "Normal 1st")
    ##lines!(ax, xplt, (1-λ) * pdf(brute_fit_dist.N₂, xplt), label = "Normal 2nd")
    ##lines!(ax, xplt, hist_curve_fit_distrib, label = "Curve fit on histogram, LsqFit.jl")
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true
    ##lines!(ax, xplt, pdf.(mixture_model_test, xplt))
    axislegend(ax)
    f
end

# ╔═╡ 377aaf8f-b909-4c42-bc77-912fd300c300
const normalization = :pdf;

# ╔═╡ 50b1a87f-49ff-4d93-aa6e-f042a87b875e
const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ 4051e244-4c84-4983-8cb9-bc7f53daa9f6
let df = CR_p_gdf_momentum[proton_momentum_index]
    f = Figure()
    ax = Axis(
        f[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Histogram of protons dN/dp at p = 10^$log_p_nat_at_slice_p mₚc",
        axis_properties...)

    if do_plot_pf
        N = df.log_dNdp_cr_pf |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "plasma frame"; bins, normalization, color = color_pf_p)
    end

    if do_plot_sf
        N = df.log_dNdp_cr_sf |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "shock frame"; bins, normalization, color = color_sf_p)
    end
    if do_plot_ISM
        N = df.log_dNdp_cr_ISM |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "ISM frame"; bins, normalization, color = color_ISM_p)
    end

    try
        axislegend(ax)
    catch e
        # axislegend has no plots to work with, because the current index doesn't have any samples. stop it complaining.
    end
    f
end

# ╔═╡ 88822f52-aab8-4931-9091-1909da6c604b
let df = CR_e_gdf_momentum[electron_momentum_index]
    f = Figure()
    ax = Axis(
        f[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Histogram of electrons dN/dp at p = 10^$log_p_nat_at_slice_e mₚc",
        axis_properties...)

    if do_plot_pf
        N = df.log_dNdp_cr_pf |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "plasma frame"; bins, normalization, color = color_pf_e)
    end

    if do_plot_sf
        N = df.log_dNdp_cr_sf |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "shock frame"; bins, normalization, color = color_sf_e)
    end
    if do_plot_ISM
        N = df.log_dNdp_cr_ISM |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "ISM frame"; bins, normalization, color = color_ISM_e)
    end

    try
        axislegend(ax)
    catch e
        # axislegend has no plots to work with, because the current index doesn't have any samples. stop it complaining.
    end
    f
end

# ╔═╡ Cell order:
# ╟─4a0e2184-0950-4b19-9b8b-061150d17ec5
# ╟─a5526239-2f05-4618-8868-0f552855d574
# ╠═f1ee2cb0-8274-11ef-0826-f55183647219
# ╟─cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
# ╠═7899ae97-fbc2-43e5-ac77-c6d725f0371e
# ╠═b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
# ╠═325d0adb-9001-4f1c-98fb-cfce042a09ca
# ╠═1e0808e1-a106-4f0e-8649-13989b8ca855
# ╠═547aad6f-32db-405d-9886-a727f1591101
# ╠═7a050dc5-7772-4933-959f-bf4fb478fc7d
# ╠═40efcd80-db38-4db3-a193-6e65ee5c4367
# ╠═3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
# ╠═fe2b3846-c753-4685-8704-e6fb50624989
# ╟─334b4ffc-1c5d-4743-88fb-ab383a3e6f80
# ╠═b544df91-fe2d-4396-892c-7faea2edd141
# ╟─4415022a-54dc-4f3d-a651-f66ae63dd051
# ╟─8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
# ╠═c159f801-b129-4919-85ef-29eedf977f14
# ╠═bdb9591b-b7ac-47e6-98bc-f18921bb64f9
# ╠═3777306e-eb41-413b-80a9-72cdc0228a94
# ╟─628130bf-da25-4799-8e5e-3d2db15b1e49
# ╟─bfc6a515-8189-487b-be08-746d865a78ae
# ╟─c5947192-0fa5-4063-8af2-74febf514b8b
# ╟─ea647872-9dc3-4fb9-9499-e396127703b2
# ╟─60ee4f38-e85f-4a2d-b17c-579531588058
# ╟─35710ad9-f2e4-487b-be19-c29500633726
# ╟─7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
# ╟─4051e244-4c84-4983-8cb9-bc7f53daa9f6
# ╟─88822f52-aab8-4931-9091-1909da6c604b
# ╟─b7a96870-784e-4ce0-830d-d245fc16e5f4
# ╟─4ac1798d-ec27-4571-9b2a-44cb432ef0d6
# ╠═f95a0d36-5dd8-4190-98c6-06e8be2ad840
# ╠═b88ef78f-6d6f-4b38-a9af-6da4f540f8c3
# ╟─f3132403-113d-4b30-9fd0-379d28ade3c7
# ╠═75e49b40-bff0-48f5-ab57-28b185f63cc9
# ╟─90e850f3-7b48-441a-92ab-1e1f6bf04e9a
# ╠═18ae83a7-98e2-4ef0-b21c-cac428146188
# ╠═e0cc631f-28a1-42db-84fc-9e7dcc9387bf
# ╠═b63ff630-624b-4e9b-bc03-dc32fd691b05
# ╠═96b36184-f98e-4b31-a2d5-1754bb40d84a
# ╠═cf870504-0f29-4354-9a4a-76971459aeba
# ╠═3596bac9-5797-40c6-a4da-cdcc1cc9a451
# ╠═7e776afd-fb12-4c90-a489-966541540599
# ╠═ac9a9859-cf1f-4084-b784-47315c0e18c1
# ╟─e3dab56a-5560-4ea8-84d7-53ce88cedc1c
# ╠═d21499d1-f010-444a-96c2-1dee378496e7
# ╠═0b1c1d4f-6ffd-423b-bf9b-31b229488038
# ╟─247d55e1-c7a4-4ccc-bce3-694f4e46dc14
# ╠═860ad43f-8683-481e-b0d2-06194ebc1af9
# ╟─7d2a121d-7e87-471b-8cc0-034796151b84
# ╠═9bc9c0f1-c3be-4944-9fec-575c4fda3ce5
# ╠═4eeb0f0b-311d-410b-b528-cbcb6f7490a7
# ╠═4d245ac7-6329-457a-970e-8a8aa23775dc
# ╠═bc44add8-d20f-4e67-ae68-7af945020d55
# ╠═6eeb6453-dee0-45eb-89ec-19b7cb2d26c1
# ╠═bd7fa9ce-3049-44d3-844f-df048003bfc5
# ╠═53fbfeb5-7993-453f-bd57-2d3c409ed46e
# ╠═e6a6cbe6-d8e4-40e1-8c29-aea7703f35a9
# ╠═a97412e6-9681-4afa-8ceb-6f37f2f6dd0b
# ╠═0aab1add-5285-4da7-b4eb-d1445b96b035
# ╠═c305f828-96c5-4839-9524-6a890a5d68fa
# ╠═d72cc184-01ce-440d-90e3-6977f9b8af7e
# ╠═df94b5b9-959c-49b8-b0b0-d8c965c61a9b
# ╠═a10b1208-3349-4a66-abc6-097bd0d5acd4
# ╟─82855970-94bf-4a67-8504-f9536f03722c
# ╠═7f0c2a67-1631-47d7-9d81-f87b44eab1c4
# ╠═c8a9815e-ed1b-44de-8d8b-9aee518cfe4e
# ╠═00f44bb7-7074-457c-aef8-566da755d748
# ╟─c1b47c6c-66ea-4014-9c04-8aa142449178
# ╠═8999b23b-4357-4655-baa6-273b218006b7
# ╠═de7bbc47-64c6-4e3c-bd2e-71f189225d52
# ╠═dbcc47ad-952e-475b-9657-f7fd280de743
# ╠═d6516ed8-0a21-4509-a1b4-34f6521ab222
# ╠═f297f99f-5f0a-4b7b-9a8b-3d7cbf37e102
# ╠═43a719e7-97d5-4e36-ba6f-3a4ec4b02463
# ╠═b822d103-e6c1-4b76-86ea-84eb84736133
# ╟─5252c2c6-969d-45c1-839c-32db557aa4b8
# ╠═5b7baac9-e657-4666-be4d-62233362aa09
# ╠═54c75433-c9e3-4d13-863e-3a1aa51f5e3e
# ╠═78e6146f-74f1-4033-9a18-1e5f8cfdd8cd
# ╠═7140dc51-55ca-437a-a6da-8812ffc35332
# ╠═8572d3a4-405c-438c-9dfc-0d37222eee9b
# ╠═2f5f0461-0101-41fb-b785-4cd96d455476
# ╠═04913ba6-6fa8-4c68-aeb3-5d9f189894ee
# ╠═5445b481-0ea6-4543-b554-7ccd257cbac6
# ╠═66a40e6b-a7f6-4b78-b5e2-f343fb02f4fe
# ╠═77733d10-ba25-4de5-b8a3-c52a4130227c
# ╠═3b6808a5-d4fc-4347-91ec-3cd389b13534
# ╠═fae42fad-f89b-4527-8bc2-8747951ea405
# ╠═ff76d92a-c68b-450e-b742-9b1706c5d310
# ╟─b9892967-520c-40e6-8cf1-3b1eb081ce04
# ╟─a08fe436-01cf-498e-8976-6e2c3173ca11
# ╠═f12b556e-5b3b-47f8-a7e6-4547e6c13d39
# ╠═788836a8-e168-4eed-b5cd-3522e43b80a6
# ╠═1360285b-8b6a-4d1c-bbb5-c6acfeddb8b6
# ╠═4153a601-06c3-4126-ace6-d354064e03f5
# ╟─8d03de5e-d344-4efd-b9af-dd5391028780
# ╠═6d5eb940-6739-4781-9dda-7433cae3cf50
# ╠═71404de8-f8b2-4d26-b7d7-41064cae1447
# ╠═cef8f0a4-0967-4e86-bfde-7fa84c474e31
# ╠═2f36a2d9-d9f5-4bf6-8fa4-c1a07532e8bb
# ╠═baf25a3a-0d13-409e-b5d9-5a1171da28b2
# ╠═4e465b9b-b2a1-42c0-ab78-ea4f620dbe30
# ╠═c21810ac-c7d7-4faf-8b2d-8985adb268da
# ╠═59a22149-3397-4e97-9f7b-5d502aacf293
# ╠═f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
# ╠═f86707a1-9d79-4df8-8798-3f7ea1d1797c
# ╠═377aaf8f-b909-4c42-bc77-912fd300c300
# ╠═50b1a87f-49ff-4d93-aa6e-f042a87b875e
