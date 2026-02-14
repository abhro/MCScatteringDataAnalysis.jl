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

# ╔═╡ d609268f-3c94-4244-9b45-8f57a21ea97d
using Revise

# ╔═╡ 7899ae97-fbc2-43e5-ac77-c6d725f0371e
using JLD2, DataFrames

# ╔═╡ b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
using Distributions

# ╔═╡ 7a050dc5-7772-4933-959f-bf4fb478fc7d
using PlutoUI

# ╔═╡ acbe8855-c586-46e3-a72b-a556df77b547
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

# ╔═╡ cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
md"""
## Import packages
"""

# ╔═╡ c3cedbde-37a4-473b-87e4-d60295362dba
md"""
## Configure notebook appearance
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
# Read data file
"""

# ╔═╡ 3bc899e5-ce26-4384-bf7b-f0bb3820f08d
CR_p_gdf_iter = load_object(datadir("Lorentz-5-processed", "dNdp-CR-protons-iteration-split.jld2"));

# ╔═╡ 80b16c96-b0f3-42a8-8544-7fbd9c06a1d9
CR_e_gdf_iter = load_object(datadir("Lorentz-5-processed", "dNdp-CR-electrons-iteration-split.jld2"));

# ╔═╡ 694a00a3-d2f5-49a3-b5d5-47fe5305329d
md"""
Enable the following cells to inspect the `GroupedDataFrame`s:
"""

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
# Plot Cosmic Ray data
"""

# ╔═╡ 3e59e047-6f47-40de-9881-f748d7f356e9
md"""
## Plotting configuration
"""

# ╔═╡ f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
axis_properties = (
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorticksvisible = true,
    yminorticksvisible = true,
)

# ╔═╡ 9aa9dffc-cb2f-4f4a-8eca-0e6046e63c11
md"""
Set a bunch of legend properties to get AlgebraOfGraphics to behave.
"""

# ╔═╡ d36bd46d-6bcd-4d97-9448-23393109e806
# options for getting legend by AlgebraOfGraphics to cooperate
legend_properties = (
    valign = :top,
    halign = :right,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    framevisible = false,
)

# ╔═╡ 50b1a87f-49ff-4d93-aa6e-f042a87b875e
color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ 3cc54622-c4e4-4c59-8828-4aa899a51e51
markersize = 5;

# ╔═╡ 113aaa26-2df1-412a-a0e1-0a5fd25e92dd
md"""
Create a method for multiplication that makes the AlgebraOfGraphics layers look a little more algebraic.
"""

# ╔═╡ 25204be4-5005-44a6-987d-f28789485a60
# A little type-piracy makes the world go round
Base.:*(b::Bool, l::Layer) = b ? l : zerolayer()

# ╔═╡ 9ad1edfe-818d-4a26-8658-632369a90845
x_map = :log_p_nat => "log p (nat)";

# ╔═╡ 968fd2bf-172b-462e-8c45-4ab7cf21f41e
visual_layer = visual(Lines);

# ╔═╡ cc0e1155-60cd-4b82-910a-9a80b489f58c
md"""
Define controllers for which iterations to plot and which frames to plot withni them
"""

# ╔═╡ 1abc4941-d902-450d-9694-aee830e6301d
plot_pf_binder = @bind do_plot_pf  CheckBox(default = true);

# ╔═╡ de79813e-1d8d-45c0-b291-203b7f19f23e
plot_sf_binder = @bind do_plot_sf  CheckBox(default = false);

# ╔═╡ c02fd290-c706-45a4-b963-c184d3bd6b2f
plot_ISM_binder = @bind do_plot_ISM CheckBox(default = false);

# ╔═╡ 3fccf366-bf6d-4c7a-a3d1-916b8f13afd3
map_layer = let
    y_label = "log(dN/dp)"

    pf_map = mapping(x_map, :log_dNdp_cr_pf => y_label, color = direct("plasma frame"))
    sf_map = mapping(x_map, :log_dNdp_cr_sf => y_label, color = direct("shock frame"))
    ISM_map = mapping(x_map, :log_dNdp_cr_ISM => y_label, color = direct("ISM frame"))

    do_plot_pf * pf_map + do_plot_sf * sf_map + do_plot_ISM * ISM_map
end;

# ╔═╡ d9b28dbe-b3d6-47d6-91c9-21b9350d5069
idx_CR_p_gdf = axes(CR_p_gdf_iter, 1);

# ╔═╡ f4be57bc-395d-4237-950d-c6d0d2b3e12c
index_binder = @bind plot_iter NumberField(idx_CR_p_gdf, default = 1);

# ╔═╡ 19a41e11-d031-498c-adbb-082e682fb67e
md"""
## Individual iterations
"""

# ╔═╡ ecf80697-b786-4b02-9563-f3d082383b76
md"""
Choose which frames to plot:
- Plasma frame: $plot_pf_binder
- Shock frame:  $plot_sf_binder
- ISM frame:    $plot_ISM_binder
"""

# ╔═╡ 338974d9-8168-4d9e-9c1e-4492bff1cf30
md"""
Select which iteration to plot:
`plot_iter` = $(index_binder) (min: $(minimum(idx_CR_p_gdf)), max: $(maximum(idx_CR_p_gdf)))
"""

# ╔═╡ bc0dcf3d-94d4-4c75-8698-13c7fb708314
CR_p_gdf_iter[plot_iter]

# ╔═╡ c2b3d96a-216e-4abe-8b0f-625419ac072f
CR_e_gdf_iter[plot_iter]

# ╔═╡ 47a47a1d-4247-4a1a-a629-0a580253b41d
md"""
at `plot_iter = 1`, there's a weird kick at the end. why?

Dr. Warren suggestion: momentum splitting. to be investigated
"""

# ╔═╡ d7d554cf-2f16-49e1-849d-25b5088e85ff
md"""
Select which iteration to plot:
`plot_iter` = $(index_binder) (min: $(minimum(idx_CR_p_gdf)), max: $(maximum(idx_CR_p_gdf)))
"""

# ╔═╡ 220c3ca5-e0b5-4f5c-86b0-e5d7cdd67558
# The primary plot in this notebook, i.e., one spectra, is through
# AlgebraOfGraphics, which is built on top of Makie.
let fig = Figure(), df = CR_p_gdf_iter[plot_iter]
    # Create a plotting specification, which specifies our data source, how to
    # transform our data (`map_layer`), and how to present the data (`visual_layer`).
    spec = data(df) * map_layer * visual_layer
    title = "dN/dp of Cosmic rays (protons), iteration $plot_iter"
    plt = draw!(fig[1, 1], spec, axis = (; title, axis_properties...))
    legend!(fig[1, 1], plt; legend_properties...)
    fig
end

# ╔═╡ 7879a41a-a284-452b-9505-a239209f1ed0
let fig = Figure(), df = CR_e_gdf_iter[plot_iter]
    spec = data(df) * map_layer * visual_layer
    title = "dN/dp of Cosmic rays (electrons), iteration $plot_iter"
    plt = draw!(fig[1, 1], spec, axis = (; title, axis_properties...))
    legend!(fig[1, 1], plt; legend_properties...)
    fig
end

# ╔═╡ 1529a53f-a084-40fc-80b0-3f9f31a5868e
md"""
Plot of ``\log(p^σ dN/dp)`` vs. ``\log(p)``
"""

# ╔═╡ 6537effb-12e6-4f4e-b34f-15dd33547921
σ = 2.23;

# ╔═╡ 79429dd7-acb7-4c5f-8843-e93e8c4ea68d
flatten_log_dNdp = (log_p, log_dNdp) -> log_p * σ + log_dNdp;

# ╔═╡ dc6d55ec-e30d-4979-9d2c-a3b6b529cdf8
md"""
Create an AlgebraOfGraphics layer which transforms the ``\log(dN/dp)`` column to ``σ \log(p) + \log(dN/dp)``.
"""

# ╔═╡ c8f1dec5-0566-4ce2-8361-21b7d1fe7480
σ_map_layer = let
    y_label = "log(dN/dp) + σ log(p)"

    pf_map = mapping(x_map, (:log_p_nat, :log_dNdp_cr_pf) => flatten_log_dNdp => y_label, color = direct("plasma frame"))
    sf_map = mapping(x_map, (:log_p_nat, :log_dNdp_cr_sf) => flatten_log_dNdp => y_label, color = direct("shock frame"))
    ISM_map = mapping(x_map, (:log_p_nat, :log_dNdp_cr_ISM) => flatten_log_dNdp => y_label, color = direct("ISM frame"))

    do_plot_pf * pf_map + do_plot_sf * sf_map + do_plot_ISM * ISM_map
end

# ╔═╡ b352849e-eca0-4ac5-acbe-9f48d0507f38
md"""
Select which iteration to plot:
`plot_iter` = $(index_binder) (min: $(minimum(idx_CR_p_gdf)), max: $(maximum(idx_CR_p_gdf)))
"""

# ╔═╡ 822bf68e-5335-45b8-b313-1cc92b53ea01
md"""
Choose which frames to plot:
- Plasma frame: $plot_pf_binder
- Shock frame:  $plot_sf_binder
- ISM frame:    $plot_ISM_binder
"""

# ╔═╡ 64eb4739-9a79-4c8c-99a3-bd338b3af6a0
let fig = Figure()
    spec = data(CR_p_gdf_iter[plot_iter]) * σ_map_layer * visual_layer
    title = "dN/dp of Cosmic rays (protons), iteration $plot_iter"
    plt = draw!(fig[1, 1], spec; axis = (; title, axis_properties...))
    legend!(fig[1, 1], plt; legend_properties..., halign = :center, valign = :bottom)

    ax = only(plt).axis # this is stupid
    hlines!(ax, 57.8, linewidth = 0.5)

    fig
end

# ╔═╡ 63db1015-02a2-4623-aa6a-b6bd772024fa
let fig = Figure()
    spec = data(CR_e_gdf_iter[plot_iter]) * σ_map_layer * visual_layer
    title = "dN/dp of Cosmic rays (electrons), iteration $plot_iter"
    plt = draw!(fig[1, 1], spec; axis = (; title, axis_properties...))
    legend!(fig[1, 1], plt; legend_properties..., halign = :center, valign = :bottom)

    ax = only(plt).axis # this is stupid
    hlines!(ax, 56.5, linewidth = 0.5)

    fig
end

# ╔═╡ 67f27108-eb9d-49b0-95ae-e016973e02b5
md"""
## Multiple iterations
"""

# ╔═╡ 6b5ff185-8eea-4e62-9cfd-3395de039b35
proton_iterations = 5620:5630;

# ╔═╡ 6c07e039-2575-49a6-a50d-531c40ee7965
let fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "dN/dp of Cosmic rays (protons)",
        xlabel = "log(p) (nat)", ylabel = "log(dN/dp) + σ log(p)",
        axis_properties...
    )

    for (i, dfp) in enumerate(CR_p_gdf_iter[proton_iterations])
        log_p, log_dNdp = dfp.log_p_nat, dfp.log_dNdp_cr_pf
        scatterlines!(ax, log_p, log_dNdp + σ * log_p, label = "plasma frame (iter $i)"; markersize)
    end

    xlims!(ax, 2, 5)
    ylims!(ax, 57.2, 58.3)
    fig
end

# ╔═╡ d76b122b-9881-4b83-ab07-34ffe17d72c3
electron_iterations = 5775:5779;

# ╔═╡ 1f35f220-7739-4097-b51d-0ab6000be247
let fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "dN/dp of Cosmic rays (electrons)",
        xlabel = "log(p) (nat)", ylabel = "log(dN/dp) + σ log(p)",
        axis_properties...
    )

    for (i, dfe) in enumerate(CR_e_gdf_iter[electron_iterations])
        log_p, log_dNdp = dfe.log_p_nat, dfe.log_dNdp_cr_pf
        scatterlines!(ax, log_p, log_dNdp + σ * log_p, label = "plasma frame (iter $i)"; markersize)
    end

    xlims!(ax, -0.3, 5)
    ylims!(ax, 56, 57.0)
    fig
end

# ╔═╡ 5d03370e-5507-41c1-94a6-24c3b5c5e9c0
let fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Flux of Cosmic rays",
        xlabel = "log(p) (mₚc)", ylabel = "log(dN/dp)",
        axis_properties...
    )

    for (i, dfp) in enumerate(CR_p_gdf_iter[proton_iterations])
        log_p, log_dNdp = dfp.log_p_nat, dfp.log_dNdp_cr_pf
        label_tup = i == 1 ? (; label = "protons") : NamedTuple()
        # lines!(ax, log_p, log_dNdp + σ*log_p; label_tup..., color = color_pf_p)
        lines!(ax, log_p, log_dNdp; label_tup..., color = color_pf_p)
    end
    for (i, dfe) in enumerate(CR_e_gdf_iter[electron_iterations])
        log_p, log_dNdp = dfe.log_p_nat, dfe.log_dNdp_cr_pf
        label_tup = i == 1 ? (; label = "electrons") : NamedTuple()
        # lines!(ax, log_p, log_dNdp + σ*log_p; label_tup..., color = color_pf_e)
        lines!(ax, log_p, log_dNdp; label_tup..., color = color_pf_e)
    end

    # xlims!(ax, 2, 5)
    # ylims!(ax, 57.2, 58.3)
    axislegend(ax, framevisible = false)
    fig
end

# ╔═╡ e5dbf380-3480-4d96-881a-8c562b5fc6ab
md"""
## All at once
"""

# ╔═╡ 4e26e9ec-b4f2-46f8-bada-945c00cb4907
# ╠═╡ disabled = true
#=╠═╡
let fig = Figure()
    ax = Axis(
        fig[1,1];
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
    fig
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
let fig = Figure()
    ax = Axis(
        fig[1,1];
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
    fig
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─f08edae2-4f29-4274-b010-07cfb3826f1e
# ╟─cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
# ╠═f1ee2cb0-8274-11ef-0826-f55183647219
# ╠═e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
# ╠═d609268f-3c94-4244-9b45-8f57a21ea97d
# ╠═7899ae97-fbc2-43e5-ac77-c6d725f0371e
# ╠═b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
# ╠═7a050dc5-7772-4933-959f-bf4fb478fc7d
# ╠═acbe8855-c586-46e3-a72b-a556df77b547
# ╠═5eaf1fbb-9dc7-40e4-87b8-8a0af299f815
# ╠═3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
# ╠═fe2b3846-c753-4685-8704-e6fb50624989
# ╟─c3cedbde-37a4-473b-87e4-d60295362dba
# ╠═b544df91-fe2d-4396-892c-7faea2edd141
# ╟─4415022a-54dc-4f3d-a651-f66ae63dd051
# ╟─8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
# ╠═3bc899e5-ce26-4384-bf7b-f0bb3820f08d
# ╠═80b16c96-b0f3-42a8-8544-7fbd9c06a1d9
# ╟─694a00a3-d2f5-49a3-b5d5-47fe5305329d
# ╠═3a2e4aee-bc90-493e-84b1-79897934f16a
# ╠═fae99b13-0b14-45c0-989f-8d0f22f0e96c
# ╟─628130bf-da25-4799-8e5e-3d2db15b1e49
# ╟─3e59e047-6f47-40de-9881-f748d7f356e9
# ╟─f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
# ╟─9aa9dffc-cb2f-4f4a-8eca-0e6046e63c11
# ╟─d36bd46d-6bcd-4d97-9448-23393109e806
# ╠═50b1a87f-49ff-4d93-aa6e-f042a87b875e
# ╠═3cc54622-c4e4-4c59-8828-4aa899a51e51
# ╟─113aaa26-2df1-412a-a0e1-0a5fd25e92dd
# ╠═25204be4-5005-44a6-987d-f28789485a60
# ╠═9ad1edfe-818d-4a26-8658-632369a90845
# ╠═3fccf366-bf6d-4c7a-a3d1-916b8f13afd3
# ╠═968fd2bf-172b-462e-8c45-4ab7cf21f41e
# ╟─cc0e1155-60cd-4b82-910a-9a80b489f58c
# ╠═1abc4941-d902-450d-9694-aee830e6301d
# ╠═de79813e-1d8d-45c0-b291-203b7f19f23e
# ╠═c02fd290-c706-45a4-b963-c184d3bd6b2f
# ╠═d9b28dbe-b3d6-47d6-91c9-21b9350d5069
# ╠═f4be57bc-395d-4237-950d-c6d0d2b3e12c
# ╟─19a41e11-d031-498c-adbb-082e682fb67e
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
# ╟─338974d9-8168-4d9e-9c1e-4492bff1cf30
# ╠═bc0dcf3d-94d4-4c75-8698-13c7fb708314
# ╠═c2b3d96a-216e-4abe-8b0f-625419ac072f
# ╟─47a47a1d-4247-4a1a-a629-0a580253b41d
# ╟─d7d554cf-2f16-49e1-849d-25b5088e85ff
# ╟─220c3ca5-e0b5-4f5c-86b0-e5d7cdd67558
# ╟─7879a41a-a284-452b-9505-a239209f1ed0
# ╟─1529a53f-a084-40fc-80b0-3f9f31a5868e
# ╠═6537effb-12e6-4f4e-b34f-15dd33547921
# ╠═79429dd7-acb7-4c5f-8843-e93e8c4ea68d
# ╟─dc6d55ec-e30d-4979-9d2c-a3b6b529cdf8
# ╠═c8f1dec5-0566-4ce2-8361-21b7d1fe7480
# ╟─b352849e-eca0-4ac5-acbe-9f48d0507f38
# ╟─822bf68e-5335-45b8-b313-1cc92b53ea01
# ╟─64eb4739-9a79-4c8c-99a3-bd338b3af6a0
# ╟─63db1015-02a2-4623-aa6a-b6bd772024fa
# ╟─67f27108-eb9d-49b0-95ae-e016973e02b5
# ╠═6b5ff185-8eea-4e62-9cfd-3395de039b35
# ╟─6c07e039-2575-49a6-a50d-531c40ee7965
# ╠═d76b122b-9881-4b83-ab07-34ffe17d72c3
# ╟─1f35f220-7739-4097-b51d-0ab6000be247
# ╟─5d03370e-5507-41c1-94a6-24c3b5c5e9c0
# ╟─e5dbf380-3480-4d96-881a-8c562b5fc6ab
# ╠═4e26e9ec-b4f2-46f8-bada-945c00cb4907
# ╟─a0be5567-9256-4c03-9a96-11d4d1973347
# ╠═526cd197-9ec7-445a-9018-3163d3916e10
