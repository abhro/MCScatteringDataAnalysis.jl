### A Pluto.jl notebook ###
# v0.20.16

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

# ╔═╡ 547aad6f-32db-405d-9886-a727f1591101
begin
    using CairoMakie
    #using WGLMakie
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

# ╔═╡ 288e5ebd-55f2-47ce-8836-98ea707bad71
using NLsolve

# ╔═╡ 2b6896b4-799e-4f00-99d7-01bbbc87e5aa
using LinearAlgebra

# ╔═╡ 75e49b40-bff0-48f5-ab57-28b185f63cc9
using BiNormalDistributions

# ╔═╡ 4d538e1f-31d3-46b4-832e-0dad096089c2
using StatsBase

# ╔═╡ ee6fd4d1-d341-4c60-bf93-8130266b5d48
using NonlinearSolve

# ╔═╡ a5526239-2f05-4618-8868-0f552855d574
md"""
# Preamble
"""

# ╔═╡ cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
md"""
## Import packages
"""

# ╔═╡ dc0952b3-5443-4c99-8cc6-497897c38dea
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

# ╔═╡ 22088c63-a7ac-45c0-97db-c1fc3360fe2d
const datadir = "G:/My Drive/MC Scattering/Processed-data";

# ╔═╡ bdb9591b-b7ac-47e6-98bc-f18921bb64f9
CR_p_gdf_momentum = load_object(joinpath(datadir, "dNdp-CR-protons-momentum-split.jld2"));

# ╔═╡ 3777306e-eb41-413b-80a9-72cdc0228a94
CR_e_gdf_momentum = load_object(joinpath(datadir, "dNdp-CR-electrons-momentum-split.jld2"));

# ╔═╡ 710739bb-10f4-4a8e-abc2-b884c6b9dfef
CR_p_gdf_momentum

# ╔═╡ 67ec1e13-d315-4566-a877-2346dca07a0c
CR_e_gdf_momentum

# ╔═╡ 628130bf-da25-4799-8e5e-3d2db15b1e49
md"""
# Plot Cosmic Ray data
"""

# ╔═╡ 067dc6df-b9ef-45cf-a32f-128a5a5cfdb6
begin
    function printstats(io::IO, gdf)
        headers = "p index | log p | nrows | # pf samples| # sf samples| # ISM samples"
        splitter = "--------|-------|-------|-------------|-------------|--------------"
        println(io, headers)
        println(io, splitter)
        for (i, df) in enumerate(gdf)
            log_p = keys(gdf)[i] |> values |> first

            println(io, @sprintf("%-7i | %5.1f | %5i | %11i | %11i | %12i", i, log_p,
                                    nrow(df),
                                    count(!ismissing, df.log_dNdp_cr_pf),
                                    count(!ismissing, df.log_dNdp_cr_sf),
                                    count(!ismissing, df.log_dNdp_cr_ISM)))
        end
        println(io, splitter)
        println(io, headers)
    end
    printstats(gdf) = printstats(stdout, gdf)
end

# ╔═╡ c16d7ba5-3272-4903-93eb-5ebf06511243
with_terminal() do
    println("For protons:")
    println()
    println()
    printstats(CR_p_gdf_momentum)
end

# ╔═╡ 653d7f99-776e-495c-8200-5e9510648de3
with_terminal() do
    println("For electrons:")
    println()
    println()
    printstats(CR_e_gdf_momentum)
end

# ╔═╡ 59a22149-3397-4e97-9f7b-5d502aacf293
const markersize = 6;

# ╔═╡ f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
const axis_properties = (xminorgridvisible = true, yminorgridvisible = true, xlabel = "log(dN/dp)");

# ╔═╡ f86707a1-9d79-4df8-8798-3f7ea1d1797c
const bins = 90;

# ╔═╡ 377aaf8f-b909-4c42-bc77-912fd300c300
const normalization = :pdf;

# ╔═╡ 50b1a87f-49ff-4d93-aa6e-f042a87b875e
const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ 35710ad9-f2e4-487b-be19-c29500633726
let
    idx_range = axes(CR_p_gdf_momentum,1)
    binder = @bind proton_momentum_index NumberField(idx_range, default = 13)
    min_idx, max_idx = extrema(idx_range)
    md"""
    Proton momentum slice to plot (index): $binder (min: $min_idx, max: $max_idx)
    """ # should the proton_momentum_index variable be considered a leak here?
end

# ╔═╡ 71404de8-f8b2-4d26-b7d7-41064cae1447
log_p_nat_at_slice = keys(CR_p_gdf_momentum)[proton_momentum_index] |> values |> first;

# ╔═╡ 7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
let
    idx_range = axes(CR_e_gdf_momentum,1)
    binder = @bind electron_momentum_index NumberField(idx_range, default = 13)
    min_idx, max_idx = extrema(idx_range)
    md"""
    Electron momentum slice to plot (index): $binder (min: $min_idx, max: $max_idx)
    """ # should the electron_momentum_index variable be considered a leak here?
end

# ╔═╡ cef8f0a4-0967-4e86-bfde-7fa84c474e31
log_p_nat_at_slice_e = keys(CR_p_gdf_momentum)[electron_momentum_index] |> values |> first;

# ╔═╡ ecf80697-b786-4b02-9563-f3d082383b76
md"""
Choose which frames to plot:
- Plasma frame: $(@bind do_plot_pf CheckBox(default=true))
- Shock frame: $(@bind do_plot_sf CheckBox(default=false))
- ISM frame: $(@bind do_plot_ISM CheckBox(default=false))
"""

# ╔═╡ 89bcb29b-0b1c-4e3a-91cb-282c05df2bc5
md"""
Value of proton momentum at slice: 10^$(log_p_nat_at_slice) _m_ₚ_c_
"""

# ╔═╡ ddc674f6-42c4-434c-8afb-c419f9752f4e
md"""
Value of electron momentum at slice: 10^$(log_p_nat_at_slice_e) *m*ₚ*c*
"""

# ╔═╡ 6d5eb940-6739-4781-9dda-7433cae3cf50
Base.:*(x::Bool, l::AoG.Layer) = x ? l : AoG.zerolayer()

# ╔═╡ 4051e244-4c84-4983-8cb9-bc7f53daa9f6
let df = CR_p_gdf_momentum[proton_momentum_index]
    f = Figure()
    ax = Axis(
        f[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Histogram of protons dN/dp at p = 10^$log_p_nat_at_slice mₚc",
        axis_properties...)

    if do_plot_pf
        N = df.log_dNdp_cr_pf |> skipmissing |> collect
        !isempty(N) && stephist!(ax, N, label = "plasma frame"; bins, normalization, color = color_pf_p)
    end

    if do_plot_sf
        N = df.log_dNdp_cr_sf |> skipmissing |> collect
        !isempty(N) && hist!(ax, N, label = "shock frame"; bins, normalization, color = color_sf_p)
    end
    if do_plot_ISM
        N = df.log_dNdp_cr_ISM |> skipmissing |> collect
        !isempty(N) && hist!(ax, N, label = "ISM frame"; bins, normalization, color = color_ISM_p)
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
        title = "Histogram of electrons dN/dp at p = 10^$log_p_nat_at_slice mₚc",
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
# Bi-normal distribution inference (method of moments)
"""

# ╔═╡ f3212b13-682f-4be4-865b-fd0f1b450aa4
Makie.update_theme!(colormap = Makie.wong_colors())

# ╔═╡ 90e850f3-7b48-441a-92ab-1e1f6bf04e9a
length(CR_p_gdf_momentum)

# ╔═╡ 18ae83a7-98e2-4ef0-b21c-cac428146188
testset_index = 67;

# ╔═╡ e0cc631f-28a1-42db-84fc-9e7dcc9387bf
CR_p_gdf_momentum[testset_index]

# ╔═╡ b63ff630-624b-4e9b-bc03-dc32fd691b05
testset = CR_p_gdf_momentum[testset_index].log_dNdp_cr_pf |> skipmissing |> collect

# ╔═╡ 96b36184-f98e-4b31-a2d5-1754bb40d84a
filter(:log_dNdp_cr_pf => >(38.7), CR_p_gdf_momentum[testset_index])

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

# ╔═╡ 3596bac9-5797-40c6-a4da-cdcc1cc9a451
testset

# ╔═╡ 6cb6018c-ca4b-43ac-b103-bb73911bd130
moms = moment.(Ref(testset), 1:8, 0)

# ╔═╡ ba79ac1a-6d27-438f-be0d-7c19c3615bae
cmoms = moment.(Ref(testset), 1:8)

# ╔═╡ 851747e7-6c04-4fd3-8ea3-4f0a742f96ec
starter = [0.5, moms[1], cmoms[2], moms[1] * 1.2, cmoms[2]]

# ╔═╡ d4e165b4-c28b-4990-8c37-4638c208abb0
Distributions.fit_mle(BiNormal, testset)

# ╔═╡ 01e2dcd7-68b4-4a65-a2c7-7e7d7c1b5b1b
# ╠═╡ disabled = true
#=╠═╡
mixture_model_test = fit_mom(BiNormal{eltype(testset)}, testset, solver = Optim.NelderMead())
  ╠═╡ =#

# ╔═╡ 951d183e-7f63-4ef2-b806-bf2c6fd94a2c
md"""
## Use NLsolve
"""

# ╔═╡ 4049368c-1c73-4a17-8ed6-aaa02e072976
md"""
## Use SciML NonlinearSolve
"""

# ╔═╡ 13b78287-4a5f-4b9a-be43-92d157373769
p = moment.(Ref(testset), 1:4, 0)

# ╔═╡ e3a95e0b-a2e4-460a-a6f9-9135c3ec35d9
md"""
Problem function for SciML specification:
"""

# ╔═╡ 77a0ffd4-459e-4cdb-8f7e-8c2bb37651bb
u₀ = zeros(5)

# ╔═╡ 529bc217-2f91-46aa-9b59-b49385b478a0
function Makie.convert_arguments(P::Type{<:AbstractPlot}, x::AbstractVector, dist::BiNormal)
    default_ptype = isdiscrete(dist) ? ScatterLines : Lines
    ptype = plottype(P, default_ptype)
    to_plotspec(ptype, convert_arguments(ptype, x, x -> (pdf(dist, x), BiNormalDistribution.componentpdfs(dist, x)...)))
end

# ╔═╡ 8d03de5e-d344-4efd-b9af-dd5391028780
md"""
# Constants and functions
"""

# ╔═╡ e780481f-ffde-407f-8dff-bc289e0ceb40
function fitdistribution(DT::Type{<:Distribution}, x::AbstractVector{Union{Missing,T}}) where {T}
    x = collect(skipmissing(x))
    if isempty(x) # don't fit to a dataset with only missings
        return missing
    end

    return Distributions.fit(DT{T}, x)
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

# ╔═╡ ab063295-d2b0-47a2-8f99-3e894f4cd646
"""
Fit `BiNormal` distribution using method of moments.
"""
function fit_mom(::Type{BiNormal{T}}, x::AbstractVector{T}; solver = nothing) where {T}

    # calculate the first 8 moments  of the given dataset
    f, s, p, q, r, u, v, w = moments(x, 8)
    @info("Computed moments", f, s, p, q, r, u, v, w)


    # initial guess should be
    # λ = 1 (all weight in first model)
    # μ₁ = <x> (mean of dataset)
    # σ₁ = √Var(x) (std. dev. of dataset)
    # μ₂ = 0?
    # σ₂ = 0?
    #θ₀ = SVector(1, f, sqrt(s - f^2), 0, 0)
    θ₀ = [1, f, sqrt(s - f^2), 0, 0]

    #prob = NonlinearLeastSquaresProblem(NonlinearFunction(fitter), θ₀, (f, s, p, q, r, u, v, w))
    prob = OptimizationProblem(OptimizationFunction(fitter, AutoZygote()), θ₀, (f, s, p, q, r, u, v, w))
    @info("Created non-linear problem", prob)

    if isnothing(solver)
        @info("Using default alg")
        sol = solve(prob)
    else
        @info("Asked to use", solver)
        sol = solve(prob, solver)
    end
    #sol = isnothing(solver) ? solve(prob) : solve(prob, solver)
    @info("Obtained solution", sol)

    return BiNormal(sol.u...)
end


# ╔═╡ 272dded6-134f-453f-8401-5f4d64379fba
"""
    mom_residuals(θ, params)

Turn the system of equations equating moments of distribution to moments of data into a function whose roots can be found.

### Arguments
- `θ`: parameters of the BiNormal distribution
- `params`: moments of the dataset
"""
function mom_residuals(θ, params)
    #(β, μ₁, rootσ₁, μ₂, rootσ₂) = θ
    #λ = constrained(β)
    #d = BiNormal(λ, μ₁, rootσ₁^2, μ₂, rootσ₂^2)
    d = BiNormal(θ...)
    (f, s, p, q, r, u, v, w) = params

    # list of residuals
    r1 = mean(d) - f
    r2 = moment(d, Val(2)) - s
    r3 = moment(d, Val(3)) - p
    r4 = moment(d, Val(4)) - q
    r5 = moment(d, Val(5)) - r
    r6 = moment(d, Val(6)) - u
    r7 = moment(d, Val(7)) - v
    r8 = moment(d, Val(8)) - w

    err = [r1, r2, r3, r4, r5, r6, r7, r8]

    return norm(err)^2
end

# ╔═╡ cfc972de-2d9c-47e4-8be2-b85d3ba44e82
s = nlsolve(x -> mom_residuals(x, moms), starter)

# ╔═╡ 9ab84473-9bfd-419a-98c0-f50d8dca9b2b
mixture_model_test = BiNormal(abs.(s.zero)...)

# ╔═╡ 1699a863-4b3d-4666-ba5e-d8d4219f5102
let f = Figure()
    ax = Axis(f[1,1], xminorgridvisible = true, yminorgridvisible = true)
    #stephist!(ax, testset; bins, normalization = :pdf, label = "stephist")
    #plot!(ax, testset_kde, label = "KDE")
    xplt = range(extrema(testset)..., length = 1000)
    lines!(ax, xplt, mixture_model_test, label = "MoM fit from nlsolve")
    #lines!(ax, xplt, pdf.(mixture_model_test, xplt))
    axislegend(ax)
    f
end

# ╔═╡ 5a01c288-0496-4228-a73b-e61d78c2e390
mean(mixture_model_test)

# ╔═╡ 02c6c3d8-0fe1-4443-b5eb-0f5275e07169
quantile(mixture_model_test, 0.3)

# ╔═╡ 5bcd406e-8851-45db-b303-5044c6dbc1d7
s.zero

# ╔═╡ 5917ab60-d980-4ac4-bdd1-796e0bcaf043
"""
    cmom_residuals(θ, params)

Turn the system of equations equating central moments of distribution to central moments of data into a function whose roots can be found.

### Arguments
- `θ`: parameters of the BiNormal distribution
- `params`: central moments of the dataset (for first moment, use the actual mean, not the first central moment, which is 0)
"""
function cmom_residuals(θ, params)
    (λ, μ₁, σ₁, μ₂, σ₂) = θ
    (f, s, p, q, r, u, v, w) = params

    error("TODO")
end

# ╔═╡ d4b39d62-1460-430d-ad8c-2da375fce2eb
"""See the docstring for the BiNormal type (defines `β(λ)`)"""
function unconstrained(λ)
    1//2 ≤ λ ≤ 1 || throw(DomainError(λ, "argument must be within [1/2,1]"))
    return log((2λ - 1)/(2 - 2λ))
end

# ╔═╡ 2df3f7c7-f196-4dae-a668-12c7116e91cb
"""See the docstring for the `BiNormal` type (defines `λ(β)`)"""
constrained(β::Real) = (2 + exp(-β))/(2 + 2exp(-β))

# ╔═╡ d8ed7e8e-2dda-4eba-95ef-e7d7a772d20c
constrained((β, μ₁, logσ₁, μ₂, logσ₂)) = (constrained(β), μ₁, exp(logσ₁), μ₂, exp(logσ₂))

# ╔═╡ 2772b072-1f67-4df6-8fc4-1243e2b79a49
function probfunc(u, p)
    # u is the parameters for the binormal distribution. sorta
    β, μ₁, logσ₁, μ₂, logσ₂ = u
    λ = constrained(β)
    d = BiNormal(λ, μ₁, exp(logσ₁), μ₂, exp(logσ₂))
    # p is the moments
    n = length(p)
    dist_moments = moment.(Ref(d), Val.(1:n))
    return dist_moments - p
end

# ╔═╡ 63cb7be2-31ee-4228-a4cc-7437af60e1fd
prob = NonlinearProblem(probfunc, u₀, p)

# ╔═╡ db901645-7047-44c1-abb6-a0195e9f6407
sciml_sol = solve(prob)

# ╔═╡ 5705ed7e-b8a3-4a16-a26d-1a439ab16b82
sciml_sol.u

# ╔═╡ 7d6d5bfd-0368-4e81-be31-c786c111bbe6
constrained(sciml_sol.u)

# ╔═╡ f34a3f2b-1e53-4ec9-b2a7-4f95c0800dc3
sciml_sol_dist = let
    β, μ₁, logσ₁, μ₂, logσ₂ = sciml_sol.u
    λ = constrained(β)
    BiNormal(λ, μ₁, exp(logσ₁), μ₂, exp(logσ₂))
end

# ╔═╡ 1d181154-ec21-47c4-8cd8-01754c05e928
plot(sciml_sol_dist)

# ╔═╡ 2a2ef8bd-a181-47bb-a02d-9a1eab0de967
BiNormalDistributions.componentpdfs(sciml_sol_dist, 3)

# ╔═╡ Cell order:
# ╠═f1ee2cb0-8274-11ef-0826-f55183647219
# ╟─a5526239-2f05-4618-8868-0f552855d574
# ╟─cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
# ╠═7899ae97-fbc2-43e5-ac77-c6d725f0371e
# ╠═b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
# ╠═547aad6f-32db-405d-9886-a727f1591101
# ╠═7a050dc5-7772-4933-959f-bf4fb478fc7d
# ╠═40efcd80-db38-4db3-a193-6e65ee5c4367
# ╠═3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
# ╠═fe2b3846-c753-4685-8704-e6fb50624989
# ╟─dc0952b3-5443-4c99-8cc6-497897c38dea
# ╠═b544df91-fe2d-4396-892c-7faea2edd141
# ╟─4415022a-54dc-4f3d-a651-f66ae63dd051
# ╟─8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
# ╠═22088c63-a7ac-45c0-97db-c1fc3360fe2d
# ╠═bdb9591b-b7ac-47e6-98bc-f18921bb64f9
# ╠═3777306e-eb41-413b-80a9-72cdc0228a94
# ╠═710739bb-10f4-4a8e-abc2-b884c6b9dfef
# ╠═67ec1e13-d315-4566-a877-2346dca07a0c
# ╟─628130bf-da25-4799-8e5e-3d2db15b1e49
# ╟─067dc6df-b9ef-45cf-a32f-128a5a5cfdb6
# ╟─c16d7ba5-3272-4903-93eb-5ebf06511243
# ╟─653d7f99-776e-495c-8200-5e9510648de3
# ╠═59a22149-3397-4e97-9f7b-5d502aacf293
# ╠═f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
# ╠═f86707a1-9d79-4df8-8798-3f7ea1d1797c
# ╠═377aaf8f-b909-4c42-bc77-912fd300c300
# ╠═50b1a87f-49ff-4d93-aa6e-f042a87b875e
# ╠═71404de8-f8b2-4d26-b7d7-41064cae1447
# ╠═cef8f0a4-0967-4e86-bfde-7fa84c474e31
# ╟─35710ad9-f2e4-487b-be19-c29500633726
# ╟─7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
# ╟─89bcb29b-0b1c-4e3a-91cb-282c05df2bc5
# ╟─ddc674f6-42c4-434c-8afb-c419f9752f4e
# ╠═6d5eb940-6739-4781-9dda-7433cae3cf50
# ╟─4051e244-4c84-4983-8cb9-bc7f53daa9f6
# ╟─88822f52-aab8-4931-9091-1909da6c604b
# ╟─b7a96870-784e-4ce0-830d-d245fc16e5f4
# ╟─4ac1798d-ec27-4571-9b2a-44cb432ef0d6
# ╠═f95a0d36-5dd8-4190-98c6-06e8be2ad840
# ╠═b88ef78f-6d6f-4b38-a9af-6da4f540f8c3
# ╟─f3132403-113d-4b30-9fd0-379d28ade3c7
# ╠═288e5ebd-55f2-47ce-8836-98ea707bad71
# ╠═2b6896b4-799e-4f00-99d7-01bbbc87e5aa
# ╠═f3212b13-682f-4be4-865b-fd0f1b450aa4
# ╠═75e49b40-bff0-48f5-ab57-28b185f63cc9
# ╠═90e850f3-7b48-441a-92ab-1e1f6bf04e9a
# ╠═18ae83a7-98e2-4ef0-b21c-cac428146188
# ╠═e0cc631f-28a1-42db-84fc-9e7dcc9387bf
# ╠═b63ff630-624b-4e9b-bc03-dc32fd691b05
# ╠═96b36184-f98e-4b31-a2d5-1754bb40d84a
# ╠═cf870504-0f29-4354-9a4a-76971459aeba
# ╠═3596bac9-5797-40c6-a4da-cdcc1cc9a451
# ╠═4d538e1f-31d3-46b4-832e-0dad096089c2
# ╠═6cb6018c-ca4b-43ac-b103-bb73911bd130
# ╠═ba79ac1a-6d27-438f-be0d-7c19c3615bae
# ╠═851747e7-6c04-4fd3-8ea3-4f0a742f96ec
# ╠═9ab84473-9bfd-419a-98c0-f50d8dca9b2b
# ╠═1699a863-4b3d-4666-ba5e-d8d4219f5102
# ╠═5a01c288-0496-4228-a73b-e61d78c2e390
# ╠═02c6c3d8-0fe1-4443-b5eb-0f5275e07169
# ╠═d4e165b4-c28b-4990-8c37-4638c208abb0
# ╠═01e2dcd7-68b4-4a65-a2c7-7e7d7c1b5b1b
# ╟─951d183e-7f63-4ef2-b806-bf2c6fd94a2c
# ╠═cfc972de-2d9c-47e4-8be2-b85d3ba44e82
# ╠═5bcd406e-8851-45db-b303-5044c6dbc1d7
# ╟─4049368c-1c73-4a17-8ed6-aaa02e072976
# ╠═ee6fd4d1-d341-4c60-bf93-8130266b5d48
# ╠═13b78287-4a5f-4b9a-be43-92d157373769
# ╟─e3a95e0b-a2e4-460a-a6f9-9135c3ec35d9
# ╠═2772b072-1f67-4df6-8fc4-1243e2b79a49
# ╠═77a0ffd4-459e-4cdb-8f7e-8c2bb37651bb
# ╠═63cb7be2-31ee-4228-a4cc-7437af60e1fd
# ╠═db901645-7047-44c1-abb6-a0195e9f6407
# ╠═5705ed7e-b8a3-4a16-a26d-1a439ab16b82
# ╠═7d6d5bfd-0368-4e81-be31-c786c111bbe6
# ╠═f34a3f2b-1e53-4ec9-b2a7-4f95c0800dc3
# ╠═1d181154-ec21-47c4-8cd8-01754c05e928
# ╠═529bc217-2f91-46aa-9b59-b49385b478a0
# ╠═2a2ef8bd-a181-47bb-a02d-9a1eab0de967
# ╟─8d03de5e-d344-4efd-b9af-dd5391028780
# ╠═84d1d644-6a5b-44eb-ab4f-3b9b7171d6fe
# ╠═e780481f-ffde-407f-8dff-bc289e0ceb40
# ╠═ab063295-d2b0-47a2-8f99-3e894f4cd646
# ╠═272dded6-134f-453f-8401-5f4d64379fba
# ╠═5917ab60-d980-4ac4-bdd1-796e0bcaf043
# ╠═d4b39d62-1460-430d-ad8c-2da375fce2eb
# ╠═2df3f7c7-f196-4dae-a668-12c7116e91cb
# ╠═d8ed7e8e-2dda-4eba-95ef-e7d7a772d20c
