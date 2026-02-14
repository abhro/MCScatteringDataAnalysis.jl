### A Pluto.jl notebook ###
# v0.20.21

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
using DrWatson

# ╔═╡ e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
@quickactivate "MCScatteringDataAnalysis"

# ╔═╡ 4c3e74c4-99d8-4d27-8787-1ea5a00e3a27
using Revise

# ╔═╡ 7899ae97-fbc2-43e5-ac77-c6d725f0371e
using JLD2, DataFrames

# ╔═╡ b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
using Distributions

# ╔═╡ 325d0adb-9001-4f1c-98fb-cfce042a09ca
using MCScatteringDataAnalysis

# ╔═╡ 1e0808e1-a106-4f0e-8649-13989b8ca855
using CairoMakie

# ╔═╡ 547aad6f-32db-405d-9886-a727f1591101
begin
    using AlgebraOfGraphics
    import AlgebraOfGraphics as AoG
end

# ╔═╡ 7a050dc5-7772-4933-959f-bf4fb478fc7d
using PlutoUI

# ╔═╡ 3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
using Missings

# ╔═╡ fe2b3846-c753-4685-8704-e6fb50624989
using Printf

# ╔═╡ 66a40e6b-a7f6-4b78-b5e2-f343fb02f4fe
using StaticArrays

# ╔═╡ 77733d10-ba25-4de5-b8a3-c52a4130227c
using LinearAlgebra

# ╔═╡ de7bbc47-64c6-4e3c-bd2e-71f189225d52
using KernelDensity

# ╔═╡ 4a0e2184-0950-4b19-9b8b-061150d17ec5
md"""
# Plot flux kernel density estimates for each momentum slice
"""

# ╔═╡ a5526239-2f05-4618-8868-0f552855d574
md"""
## Preamble
"""

# ╔═╡ 29b0ffcc-1799-491b-9853-7296c68483cf
md"""
!!! warning "Do not run"

    The following two cells with DrWatson are commented out because the root Project.toml file no longer contain some of the packages used in ths notebook for perfomance reasons. If you need to re-run, uncomment the lines invloving `DrWatson` and `@quickactivate`, and re-add the packages to the Project.toml file.
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

# ╔═╡ bdb9591b-b7ac-47e6-98bc-f18921bb64f9
CR_p_gdf_momentum = load_object(datadir("Lorentz-5-processed", "dNdp-CR-protons-momentum-split.jld2"));

# ╔═╡ 3777306e-eb41-413b-80a9-72cdc0228a94
CR_e_gdf_momentum = load_object(datadir("Lorentz-5-processed", "dNdp-CR-electrons-momentum-split.jld2"));

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

# ╔═╡ 628130bf-da25-4799-8e5e-3d2db15b1e49
md"""
## Plot Cosmic Ray data
"""

# ╔═╡ 2af8ede1-28cc-45f7-ba86-7632b4931c17
md"""
### Plotting configurations
"""

# ╔═╡ 2f36a2d9-d9f5-4bf6-8fa4-c1a07532e8bb
proton_indices = axes(CR_p_gdf_momentum, 1);

# ╔═╡ baf25a3a-0d13-409e-b5d9-5a1171da28b2
electron_indices = axes(CR_e_gdf_momentum, 1);

# ╔═╡ 4e465b9b-b2a1-42c0-ab78-ea4f620dbe30
proton_index_binder = @bind proton_momentum_index NumberField(proton_indices, default = 13);

# ╔═╡ c21810ac-c7d7-4faf-8b2d-8985adb268da
electron_index_binder = @bind electron_momentum_index NumberField(electron_indices, default = 13);

# ╔═╡ 59a22149-3397-4e97-9f7b-5d502aacf293
markersize = 6;

# ╔═╡ f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
axis_properties = (xminorgridvisible = true, yminorgridvisible = true, xlabel = "log(dN/dp)")

# ╔═╡ f86707a1-9d79-4df8-8798-3f7ea1d1797c
bins = 90;

# ╔═╡ 377aaf8f-b909-4c42-bc77-912fd300c300
normalization = :pdf;

# ╔═╡ 50b1a87f-49ff-4d93-aa6e-f042a87b875e
color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ ecf80697-b786-4b02-9563-f3d082383b76
md"""
Choose which frames to plot:
- Plasma frame: $(@bind do_plot_pf CheckBox(default=true))
- Shock frame: $(@bind do_plot_sf CheckBox(default=false))
- ISM frame: $(@bind do_plot_ISM CheckBox(default=false))
"""

# ╔═╡ c1b47c6c-66ea-4014-9c04-8aa142449178
md"""
### Kernel density estimates
"""

# ╔═╡ 7b238d92-5588-43af-a861-ef89f2e8178d
testset = CR_p_gdf_momentum[proton_momentum_index].log_dNdp_cr_pf

# ╔═╡ 8999b23b-4357-4655-baa6-273b218006b7
Makie.density(testset)

# ╔═╡ dbcc47ad-952e-475b-9657-f7fd280de743
# ╠═╡ disabled = true
#=╠═╡
let
    fig, ax, _ = plot(testset_kde, label = "KDE")
    plot!(ax, testset_kde.x[density_maxes.indices], density_maxes.heights, label = "modes")
    axislegend(ax)
    fig
end
  ╠═╡ =#

# ╔═╡ d6516ed8-0a21-4509-a1b4-34f6521ab222
# ╠═╡ disabled = true
#=╠═╡
testset_kde.x[density_maxes.indices]
  ╠═╡ =#

# ╔═╡ 5252c2c6-969d-45c1-839c-32db557aa4b8
md"""
### Fit Kernel density estimate curve
"""

# ╔═╡ 5b7baac9-e657-4666-be4d-62233362aa09
# ╠═╡ disabled = true
#=╠═╡
testset_kde.x |> length
  ╠═╡ =#

# ╔═╡ 54c75433-c9e3-4d13-863e-3a1aa51f5e3e
# ╠═╡ disabled = true
#=╠═╡
testset_kde.density |> length
  ╠═╡ =#

# ╔═╡ 8572d3a4-405c-438c-9dfc-0d37222eee9b
# ╠═╡ disabled = true
#=╠═╡
let fig = Figure()
    ax = Axis(fig[1, 1], xminorgridvisible = true, yminorgridvisible = true)
    ##stephist!(ax, testset; bins, normalization = :pdf)
    plot!(ax, testset_kde, label = "kde")
    xplt = range(extrema(testset)..., length = 1000)
    lines!(ax, xplt, kde_curve_fit_distrib, label = "Curve fit on kde, LsqFit.jl")
    #lines!(ax, xplt, pdf.(mixture_model_test, xplt))
    axislegend(ax)
    fig
end
  ╠═╡ =#

# ╔═╡ 8d03de5e-d344-4efd-b9af-dd5391028780
md"""
## Constants and functions
"""

# ╔═╡ 6d5eb940-6739-4781-9dda-7433cae3cf50
Base.:*(x::Bool, l::AoG.Layer) = x ? l : AoG.zerolayer()

# ╔═╡ 71404de8-f8b2-4d26-b7d7-41064cae1447
log_p_nat_at_slice_p = keys(CR_p_gdf_momentum)[proton_momentum_index] |> values |> only;

# ╔═╡ 35710ad9-f2e4-487b-be19-c29500633726
md"""
Proton momentum slice to plot (index): $(proton_index_binder) (min: $(minimum(proton_indices)), max: $(maximum(proton_indices)))

Value of proton momentum at slice: 10^$(log_p_nat_at_slice_p) _m_ₚ_c_
"""

# ╔═╡ cef8f0a4-0967-4e86-bfde-7fa84c474e31
log_p_nat_at_slice_e = keys(CR_p_gdf_momentum)[electron_momentum_index] |> values |> only;

# ╔═╡ 7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
md"""
Electron momentum slice to plot (index): $electron_index_binder (min: $(minimum(electron_indices)), max: $(maximum(electron_indices)))

Value of electron momentum at slice: 10^$(log_p_nat_at_slice_e) *m*ₚ*c*
"""

# ╔═╡ Cell order:
# ╟─4a0e2184-0950-4b19-9b8b-061150d17ec5
# ╟─a5526239-2f05-4618-8868-0f552855d574
# ╟─29b0ffcc-1799-491b-9853-7296c68483cf
# ╠═f1ee2cb0-8274-11ef-0826-f55183647219
# ╠═e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
# ╠═4c3e74c4-99d8-4d27-8787-1ea5a00e3a27
# ╟─cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
# ╠═7899ae97-fbc2-43e5-ac77-c6d725f0371e
# ╠═b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
# ╠═325d0adb-9001-4f1c-98fb-cfce042a09ca
# ╠═1e0808e1-a106-4f0e-8649-13989b8ca855
# ╠═547aad6f-32db-405d-9886-a727f1591101
# ╠═7a050dc5-7772-4933-959f-bf4fb478fc7d
# ╠═3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
# ╠═fe2b3846-c753-4685-8704-e6fb50624989
# ╠═66a40e6b-a7f6-4b78-b5e2-f343fb02f4fe
# ╠═77733d10-ba25-4de5-b8a3-c52a4130227c
# ╟─334b4ffc-1c5d-4743-88fb-ab383a3e6f80
# ╠═b544df91-fe2d-4396-892c-7faea2edd141
# ╟─4415022a-54dc-4f3d-a651-f66ae63dd051
# ╟─8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
# ╠═bdb9591b-b7ac-47e6-98bc-f18921bb64f9
# ╠═3777306e-eb41-413b-80a9-72cdc0228a94
# ╟─bfc6a515-8189-487b-be08-746d865a78ae
# ╟─c5947192-0fa5-4063-8af2-74febf514b8b
# ╟─ea647872-9dc3-4fb9-9499-e396127703b2
# ╟─60ee4f38-e85f-4a2d-b17c-579531588058
# ╟─628130bf-da25-4799-8e5e-3d2db15b1e49
# ╟─2af8ede1-28cc-45f7-ba86-7632b4931c17
# ╠═2f36a2d9-d9f5-4bf6-8fa4-c1a07532e8bb
# ╠═baf25a3a-0d13-409e-b5d9-5a1171da28b2
# ╠═4e465b9b-b2a1-42c0-ab78-ea4f620dbe30
# ╠═c21810ac-c7d7-4faf-8b2d-8985adb268da
# ╠═59a22149-3397-4e97-9f7b-5d502aacf293
# ╟─f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
# ╠═f86707a1-9d79-4df8-8798-3f7ea1d1797c
# ╠═377aaf8f-b909-4c42-bc77-912fd300c300
# ╠═50b1a87f-49ff-4d93-aa6e-f042a87b875e
# ╟─35710ad9-f2e4-487b-be19-c29500633726
# ╟─7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
# ╟─c1b47c6c-66ea-4014-9c04-8aa142449178
# ╠═7b238d92-5588-43af-a861-ef89f2e8178d
# ╠═8999b23b-4357-4655-baa6-273b218006b7
# ╠═de7bbc47-64c6-4e3c-bd2e-71f189225d52
# ╠═dbcc47ad-952e-475b-9657-f7fd280de743
# ╠═d6516ed8-0a21-4509-a1b4-34f6521ab222
# ╟─5252c2c6-969d-45c1-839c-32db557aa4b8
# ╠═5b7baac9-e657-4666-be4d-62233362aa09
# ╠═54c75433-c9e3-4d13-863e-3a1aa51f5e3e
# ╠═8572d3a4-405c-438c-9dfc-0d37222eee9b
# ╟─8d03de5e-d344-4efd-b9af-dd5391028780
# ╠═6d5eb940-6739-4781-9dda-7433cae3cf50
# ╠═71404de8-f8b2-4d26-b7d7-41064cae1447
# ╠═cef8f0a4-0967-4e86-bfde-7fa84c474e31
