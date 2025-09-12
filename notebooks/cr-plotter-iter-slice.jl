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
# ╠═╡ skip_as_script = true
#=╠═╡
import Pkg; Pkg.activate(Base.current_project())
  ╠═╡ =#

# ╔═╡ 7899ae97-fbc2-43e5-ac77-c6d725f0371e
using JLD2, DataFrames

# ╔═╡ b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
using Distributions

# ╔═╡ 7a050dc5-7772-4933-959f-bf4fb478fc7d
using PlutoUI

# ╔═╡ 40efcd80-db38-4db3-a193-6e65ee5c4367
using PlutoUI: Slider

# ╔═╡ acbe8855-c586-46e3-a72b-a556df77b547
#using WGLMakie
using CairoMakie

# ╔═╡ 5eaf1fbb-9dc7-40e4-87b8-8a0af299f815
begin
    using AlgebraOfGraphics
    import AlgebraOfGraphics as AoG
end

# ╔═╡ 3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
using Missings

# ╔═╡ fe2b3846-c753-4685-8704-e6fb50624989
using Printf

# ╔═╡ f08edae2-4f29-4274-b010-07cfb3826f1e
md"""
# Plot spectra for each iteration
"""

# ╔═╡ a5526239-2f05-4618-8868-0f552855d574
md"""
## Preamble
"""

# ╔═╡ cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
md"""
### Import packages
"""

# ╔═╡ c3cedbde-37a4-473b-87e4-d60295362dba
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

# ╔═╡ ee3eab6d-7913-4650-a5c3-aabf0747a58a
const datadir = "G:/My Drive/MC Scattering/Run-data/Lorentz-5-processed";

# ╔═╡ 3bc899e5-ce26-4384-bf7b-f0bb3820f08d
CR_p_gdf_iter = load_object(joinpath(datadir, "dNdp-CR-protons-iteration-split.jld2"));

# ╔═╡ 80b16c96-b0f3-42a8-8544-7fbd9c06a1d9
CR_e_gdf_iter = load_object(joinpath(datadir, "dNdp-CR-electrons-iteration-split.jld2"));

# ╔═╡ 3a2e4aee-bc90-493e-84b1-79897934f16a
# ╠═╡ disabled = true
#=╠═╡
CR_p_gdf_iter
  ╠═╡ =#

# ╔═╡ fae99b13-0b14-45c0-989f-8d0f22f0e96c
# ╠═╡ disabled = true
#=╠═╡
CR_e_gdf_iter
  ╠═╡ =#

# ╔═╡ 628130bf-da25-4799-8e5e-3d2db15b1e49
md"""
## Plot Cosmic Ray data
"""

# ╔═╡ 3e59e047-6f47-40de-9881-f748d7f356e9
md"""
### Plotting configuration
"""

# ╔═╡ f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
const axis_properties = (
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorticksvisible = true,
    yminorticksvisible = true,
)

# ╔═╡ d36bd46d-6bcd-4d97-9448-23393109e806
# options for getting legend by AlgebraOfGraphics to cooperate
const legend_properties = (
    valign = :top,
    halign = :right,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    framevisible = false,
)

# ╔═╡ 50b1a87f-49ff-4d93-aa6e-f042a87b875e
const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ 3cc54622-c4e4-4c59-8828-4aa899a51e51
const markersize = 5;

# ╔═╡ 968fd2bf-172b-462e-8c45-4ab7cf21f41e
visual_layer = visual(Lines);

# ╔═╡ ecf80697-b786-4b02-9563-f3d082383b76
md"""
Choose which frames to plot:
- Plasma frame: $(@bind do_plot_pf  CheckBox(default=true))
- Shock frame:  $(@bind do_plot_sf  CheckBox(default=false))
- ISM frame:    $(@bind do_plot_ISM CheckBox(default=false))
"""

# ╔═╡ 3fccf366-bf6d-4c7a-a3d1-916b8f13afd3
map_layer = let
    x_map = :log_p_nat => "log p (nat)"
    y_label = "log(dN/dp)"

    # A little type-piracy makes the world go round
    Base.:*(b::Bool, l::Layer) = b ? l : zerolayer()

    pf_map = mapping(x_map, :log_dNdp_cr_pf => y_label, color = direct("plasma frame"))
    sf_map = mapping(x_map, :log_dNdp_cr_sf => y_label, color = direct("shock frame"))
    ISM_map = mapping(x_map, :log_dNdp_cr_ISM => y_label, color = direct("ISM frame"))

    do_plot_pf*pf_map + do_plot_sf*sf_map + do_plot_ISM*ISM_map
end;

# ╔═╡ d9b28dbe-b3d6-47d6-91c9-21b9350d5069
const idx_CR_p_gdf = axes(CR_p_gdf_iter, 1);

# ╔═╡ f4be57bc-395d-4237-950d-c6d0d2b3e12c
const index_binder = @bind plot_iter NumberField(idx_CR_p_gdf, default = 1);

# ╔═╡ d7d554cf-2f16-49e1-849d-25b5088e85ff
md"""
Select which iteration to plot:

`plot_iter` = $(index_binder)
"""

# ╔═╡ 19a41e11-d031-498c-adbb-082e682fb67e
md"""
### Individual iterations
"""

# ╔═╡ c2b3d96a-216e-4abe-8b0f-625419ac072f
# ╠═╡ disabled = true
#=╠═╡
CR_e_gdf_iter[plot_iter]
  ╠═╡ =#

# ╔═╡ 47a47a1d-4247-4a1a-a629-0a580253b41d
md"""
at `plot_iter = 1`, there's a weird kick at the end. why?

Dr. Warren suggestion: momentum splitting. to be investigated
"""

# ╔═╡ 220c3ca5-e0b5-4f5c-86b0-e5d7cdd67558
let f = Figure(), df = CR_p_gdf_iter[plot_iter]
    spec = data(df) * map_layer * visual_layer
    title = "dN/dp of Cosmic rays (protons), iteration $plot_iter"
    plt = draw!(f[1,1], spec, axis = (; title, axis_properties...))
    legend!(f[1,1], plt; legend_properties...)
    f
end

# ╔═╡ 7879a41a-a284-452b-9505-a239209f1ed0
let f = Figure(), df = CR_e_gdf_iter[plot_iter]
    spec = data(df) * map_layer * visual_layer
    title = "dN/dp of Cosmic rays (electrons), iteration $plot_iter"
    plt = draw!(f[1,1], spec, axis = (; title, axis_properties...))
    legend!(f[1,1], plt; legend_properties...)
    f
end

# ╔═╡ b352849e-eca0-4ac5-acbe-9f48d0507f38
md"""
Select which iteration to plot:

`plot_iter` = $(index_binder)
"""

# ╔═╡ 1529a53f-a084-40fc-80b0-3f9f31a5868e
md"""
Plot of ``\log(p^σ dN/dp)`` vs. ``\log(p)``
"""

# ╔═╡ 6537effb-12e6-4f4e-b34f-15dd33547921
const σ = 2.23;

# ╔═╡ f4930314-a64c-4b6a-bcef-c0d9dcf2ef81
let
    dfp = CR_p_gdf_iter[plot_iter]
    dfe = CR_e_gdf_iter[plot_iter]

    f = Figure()

    ax = Axis(
        f[1,1]; title = "dN/dp of Cosmic rays, iteration $plot_iter",
        xlabel = "log(p) (nat)", ylabel = "log(dN/dp) + σ log(p)", axis_properties...)

    if do_plot_pf
        log_p, log_dNdp = (dfp.log_p_nat, dfp.log_dNdp_cr_pf)
        scatterlines!(ax, log_p, log_dNdp .+ σ*log_p, label = "protons, plasma frame"; color = color_pf_p, markersize)

        log_p, log_dNdp = (dfe.log_p_nat, dfe.log_dNdp_cr_pf)
        scatterlines!(ax, log_p, log_dNdp .+ σ*log_p, label = "electrons, plasma frame"; color = color_pf_e, markersize)
    end
    if do_plot_sf
        log_p, log_dNdp = (dfp.log_p_nat, dfp.log_dNdp_cr_sf)
        scatterlines!(ax, log_p, log_dNdp .+ σ*log_p, label = "protons, shock frame"; color = color_sf_p, markersize)

        log_p, log_dNdp = (dfe.log_p_nat, dfe.log_dNdp_cr_sf)
        scatterlines!(ax, log_p, log_dNdp .+ σ*log_p, label = "electrons, shock frame"; color = color_sf_e, markersize)
    end
    if do_plot_ISM
        log_p, log_dNdp = (dfp.log_p_nat, dfp.log_dNdp_cr_ISM)
        scatterlines!(ax, log_p, log_dNdp .+ σ*log_p, label = "protons, ISM frame"; color = color_ISM_p, markersize)

        log_p, log_dNdp = (dfe.log_p_nat, dfe.log_dNdp_cr_ISM)
        scatterlines!(ax, log_p, log_dNdp .+ σ*log_p, label = "electrons, ISM frame"; color = color_ISM_e, markersize)
    end

    hlines!(ax, 57.8, color = color_pf_p, linewidth = 0.5)
    hlines!(ax, 56.5, color = color_pf_e, linewidth = 0.5)

    #xlims!(ax, -1, 8)
    #ylims!(ax, 56, 58.5)
    try
        axislegend(ax, position = :cb)
    catch e
        @error(e)
    end
    f
end

# ╔═╡ 67f27108-eb9d-49b0-95ae-e016973e02b5
md"""
### Multiple iterations
"""

# ╔═╡ 6c07e039-2575-49a6-a50d-531c40ee7965
let f = Figure()
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (protons)",
        xlabel = "log(p) (nat)", ylabel = "log(dN/dp) + σ log(p)",
        axis_properties...)

    for (i, dfp) in enumerate(CR_p_gdf_iter[5620:5630])
        log_p, log_dNdp = dfp.log_p_nat, dfp.log_dNdp_cr_pf
        scatterlines!(ax, log_p, log_dNdp + σ*log_p, label = "plasma frame (iter $i)"; markersize)
    end

    xlims!(ax, 2, 5)
    ylims!(ax, 57.2, 58.3)
    f
end

# ╔═╡ 1f35f220-7739-4097-b51d-0ab6000be247
let f = Figure()
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (electrons)",
        xlabel = "log(p) (nat)", ylabel = "log(dN/dp) + σ log(p)",
        axis_properties...)

    for (i, dfe) in enumerate(CR_e_gdf_iter[5775:5779])
        log_p, log_dNdp = dfe.log_p_nat, dfe.log_dNdp_cr_pf
        scatterlines!(ax, log_p, log_dNdp + σ*log_p, label = "plasma frame (iter $i)"; markersize)
    end

    xlims!(ax, -0.3, 5)
    ylims!(ax, 56, 57.0)
    f
end

# ╔═╡ e5dbf380-3480-4d96-881a-8c562b5fc6ab
md"""
### All at once
"""

# ╔═╡ 4e26e9ec-b4f2-46f8-bada-945c00cb4907
# ╠═╡ disabled = true
#=╠═╡
let f = Figure()
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (protons)",
        xlabel = "log(p) (nat)", ylabel = "log(dN/dp) + σ log(p)",
        axis_properties...)

    for (i, dfp) in enumerate(CR_p_gdf_iter)
        log_p, log_dNdp = dfp.log_p_nat, dfp.log_dNdp_cr_pf
        scatterlines!(ax, log_p, log_dNdp + σ*log_p, label = "plasma frame (iter $i)"; markersize)
    end

    #hlines!(ax, 57.5)

    #xlims!(ax, 2, 5)
    #ylims!(ax, 57.25, 58.7)
    #axislegend(ax, position = :lb)
    f
end
  ╠═╡ =#

# ╔═╡ a0be5567-9256-4c03-9a96-11d4d1973347
md"""
Find the iteration number of the anomalous electron run

Found: iteration 5775 (iseed: 289, iter within seed: 15).
"""

# ╔═╡ 526cd197-9ec7-445a-9018-3163d3916e10
# ╠═╡ disabled = true
#=╠═╡
let f = Figure()
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (electrons)",
        xlabel = "log(p) (nat)", ylabel = "log(dN/dp) + σ log(p)",
        axis_properties...)

    for (i, dfe) in enumerate(CR_e_gdf_iter)
        log_p, log_dNdp = dfe.log_p_nat, dfe.log_dNdp_cr_pf
        scatterlines!(ax, log_p, log_dNdp + σ*log_p, label = "plasma frame (iter $i)"; markersize)
    end

    xlims!(ax, -0.3, 5)
    ylims!(ax, 56, 57.4)
    #axislegend(ax, position = :lb)
    f
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─f08edae2-4f29-4274-b010-07cfb3826f1e
# ╟─a5526239-2f05-4618-8868-0f552855d574
# ╠═f1ee2cb0-8274-11ef-0826-f55183647219
# ╟─cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
# ╠═7899ae97-fbc2-43e5-ac77-c6d725f0371e
# ╠═b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
# ╠═7a050dc5-7772-4933-959f-bf4fb478fc7d
# ╠═40efcd80-db38-4db3-a193-6e65ee5c4367
# ╠═acbe8855-c586-46e3-a72b-a556df77b547
# ╠═5eaf1fbb-9dc7-40e4-87b8-8a0af299f815
# ╠═3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
# ╠═fe2b3846-c753-4685-8704-e6fb50624989
# ╟─c3cedbde-37a4-473b-87e4-d60295362dba
# ╠═b544df91-fe2d-4396-892c-7faea2edd141
# ╟─4415022a-54dc-4f3d-a651-f66ae63dd051
# ╟─8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
# ╠═ee3eab6d-7913-4650-a5c3-aabf0747a58a
# ╠═3bc899e5-ce26-4384-bf7b-f0bb3820f08d
# ╠═80b16c96-b0f3-42a8-8544-7fbd9c06a1d9
# ╠═3a2e4aee-bc90-493e-84b1-79897934f16a
# ╠═fae99b13-0b14-45c0-989f-8d0f22f0e96c
# ╟─628130bf-da25-4799-8e5e-3d2db15b1e49
# ╟─3e59e047-6f47-40de-9881-f748d7f356e9
# ╟─f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
# ╟─d36bd46d-6bcd-4d97-9448-23393109e806
# ╠═50b1a87f-49ff-4d93-aa6e-f042a87b875e
# ╠═3cc54622-c4e4-4c59-8828-4aa899a51e51
# ╠═3fccf366-bf6d-4c7a-a3d1-916b8f13afd3
# ╠═968fd2bf-172b-462e-8c45-4ab7cf21f41e
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
# ╠═d9b28dbe-b3d6-47d6-91c9-21b9350d5069
# ╠═f4be57bc-395d-4237-950d-c6d0d2b3e12c
# ╟─19a41e11-d031-498c-adbb-082e682fb67e
# ╠═c2b3d96a-216e-4abe-8b0f-625419ac072f
# ╟─47a47a1d-4247-4a1a-a629-0a580253b41d
# ╟─d7d554cf-2f16-49e1-849d-25b5088e85ff
# ╟─220c3ca5-e0b5-4f5c-86b0-e5d7cdd67558
# ╟─7879a41a-a284-452b-9505-a239209f1ed0
# ╟─b352849e-eca0-4ac5-acbe-9f48d0507f38
# ╟─1529a53f-a084-40fc-80b0-3f9f31a5868e
# ╠═6537effb-12e6-4f4e-b34f-15dd33547921
# ╟─f4930314-a64c-4b6a-bcef-c0d9dcf2ef81
# ╟─67f27108-eb9d-49b0-95ae-e016973e02b5
# ╟─6c07e039-2575-49a6-a50d-531c40ee7965
# ╟─1f35f220-7739-4097-b51d-0ab6000be247
# ╟─e5dbf380-3480-4d96-881a-8c562b5fc6ab
# ╠═4e26e9ec-b4f2-46f8-bada-945c00cb4907
# ╟─a0be5567-9256-4c03-9a96-11d4d1973347
# ╠═526cd197-9ec7-445a-9018-3163d3916e10
