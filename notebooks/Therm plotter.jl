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

# ╔═╡ b0d4ca90-bb62-11f0-80dd-b560b36acdf6
using DrWatson

# ╔═╡ 63875dd0-d278-44c4-a85a-78ca6d7e9e10
@quickactivate "MCScatteringDataAnalysis"

# ╔═╡ a08cb144-1980-46e6-b835-65c93fd5a265
using Revise

# ╔═╡ d3e5bcd8-31e0-424f-8626-5eb4b679b818
using PlutoUI

# ╔═╡ b1833ac8-965e-45d5-bffd-8e31fac39939
using JLD2

# ╔═╡ 24564257-76f5-4200-b9a8-30180e853233
using DataFrames

# ╔═╡ e50e35b0-885a-4f34-b716-a0167f62d8f7
using CairoMakie

# ╔═╡ 9e779980-f661-45cd-bdd1-20e073df9ede
begin
    using AlgebraOfGraphics
    import AlgebraOfGraphics as AoG
end

# ╔═╡ a1e20290-55e9-4355-bc45-74c5f825565f
md"""
# Plot thermal data
"""

# ╔═╡ d6dbd935-f720-464c-a3ed-9a8f26b8d287
md"""
## Configure notebook appearance
"""

# ╔═╡ 4f092e96-b696-4a79-a109-95c2ecee448f
TableOfContents(depth = 6)

# ╔═╡ e368d071-20c3-433d-b884-c0ecb3887b5a
# Increase cell width
html"""<style>
main {
    max-width: 83%;
    padding-left: max(300px, 5%);
    padding-right: 0%;
}
</style>"""

# ╔═╡ 537c5e04-31ef-4c5e-a80e-b255ed166ba1
md"""
# Read data file
"""

# ╔═╡ 9b3120c2-2b26-41e5-a067-19cadd0c5136
therm_p_gdf_iter_filename = datadir("Lorentz-5-processed", "dNdp-therm-protons-iteration-split.jld2");

# ╔═╡ 7bc31948-1d91-461d-8fa0-3a58aebeace4
therm_e_gdf_iter_filename = datadir("Lorentz-5-processed", "dNdp-therm-electrons-iteration-split.jld2");

# ╔═╡ b6a0acf6-628c-40c5-b5f9-c47ba3814ab2
therm_p_gdf_iter = load_object(therm_p_gdf_iter_filename);

# ╔═╡ 667ae8fa-015e-4c30-b8b0-f1bc9646a5ab
therm_e_gdf_iter = load_object(therm_e_gdf_iter_filename);

# ╔═╡ 4e419982-323d-4cdf-b026-4ebec2cd0c39
md"""
Enable the following cells to inspect the `GroupedDataFrame`s:
"""

# ╔═╡ fa5460cd-c3e5-440c-8c13-94f8374fa643
# ╠═╡ disabled = true
#=╠═╡
therm_p_gdf_iter
  ╠═╡ =#

# ╔═╡ 6764b6b6-63cf-4a7b-8750-65002ee0a728
# ╠═╡ disabled = true
#=╠═╡
therm_e_gdf_iter
  ╠═╡ =#

# ╔═╡ c8442393-614f-4bfa-aeab-c8fad1b3e0b8
md"""
# Plot thermal particle data
"""

# ╔═╡ 34fb2a5f-6c38-43a7-aec5-0c05d163ab4e
md"""
## Plotting configuration
"""

# ╔═╡ 3dd36273-6ab6-4003-a7f3-4ca4a456ed21
axis_properties = (
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorticksvisible = true,
    yminorticksvisible = true,
)

# ╔═╡ 76507571-215c-457f-8c5f-affdca267b1f
md"""
Set a bunch of legend properties to get AlgebraOfGraphics to behave.
"""

# ╔═╡ 01f0eec1-3f3e-44b3-b9cd-246eb3b9e557
# options for getting legend by AlgebraOfGraphics to cooperate
legend_properties = (
    valign = :top,
    halign = :right,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    framevisible = false,
)

# ╔═╡ c14fef4c-6dcd-410c-8d23-8beee9a38e1e
color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ 10e30fc9-7a25-4672-86e9-724d490a0db4
markersize = 5;

# ╔═╡ 0ce3bbd3-150f-4fb5-a7f6-4264b11dfbe6
md"""
Create a method for multiplication that makes the AlgebraOfGraphics layers look a little more algebraic.
"""

# ╔═╡ d5fe4e4c-c089-43fb-a8cc-26053ef5775f
# A little type-piracy makes the world go round
Base.:*(b::Bool, l::Layer) = b ? l : zerolayer()

# ╔═╡ c4d125ac-7412-4f52-838c-33e1bfa23386
pf_map = mapping(:log_dNdp_therm_pvals_nat_pf => "log p (nat)", :log_dNdp_therm_pf => "log(dN/dp)", color = direct("plasma frame"));

# ╔═╡ 25858c85-8006-4572-ad81-e8b220fab0d8
sf_map = mapping(:log_dNdp_therm_pvals_nat_sf => "log p (nat)", :log_dNdp_therm_sf => "log(dN/dp)", color = direct("shock frame"));

# ╔═╡ 4e8d3b62-0c32-427b-b1f2-e1a74d82a66c
ISM_map = mapping(:log_dNdp_therm_pvals_nat_ISM => "log p (nat)", :log_dNdp_therm_ISM => "log(dN/dp)", color = direct("ISM frame"));

# ╔═╡ da63beaa-882c-4d34-b69a-3423d98e12c1
visual_layer = visual(Lines);

# ╔═╡ 78db62d2-4717-4bcf-a6fb-6cd80a2b9faf
md"""
Define controllers for which iterations to plot and which frames to plot withni them
"""

# ╔═╡ 2e5e7272-0313-471d-8987-de64f79b84b0
plot_pf_binder = @bind do_plot_pf  CheckBox(default = true);

# ╔═╡ 8c12b2c3-b5e3-4145-97b8-bba2ab5f03fc
plot_sf_binder = @bind do_plot_sf  CheckBox(default = false);

# ╔═╡ dc4d6726-25e1-4aeb-9029-fa88282667db
plot_ISM_binder = @bind do_plot_ISM CheckBox(default = false);

# ╔═╡ 3b1ecd57-8f9f-40f7-818c-651a8a320f6f
map_layer = do_plot_pf * pf_map + do_plot_sf * sf_map + do_plot_ISM * ISM_map;

# ╔═╡ 6aae8a53-9e2a-4c54-b5da-af9452fa2d5b
idx_therm_p_gdf = axes(therm_p_gdf_iter, 1);

# ╔═╡ 911db519-1f6c-4534-9de7-5e03e0869ae3
index_binder = @bind plot_iter NumberField(idx_therm_p_gdf, default = 1);

# ╔═╡ 720c2380-c92a-4a1f-a8f4-5e455a79651a
md"""
## Individual iterations
"""

# ╔═╡ 6c256957-03b3-40e6-9961-e26568c2c72c
md"""
Choose which frames to plot:
- Plasma frame: $plot_pf_binder
- Shock frame:  $plot_sf_binder
- ISM frame:    $plot_ISM_binder
"""

# ╔═╡ f3e07ca1-3e4c-41fc-b877-dc635097b15e
md"""
Select which iteration to plot:
`plot_iter` = $(index_binder) (min: $(minimum(idx_therm_p_gdf)), max: $(maximum(idx_therm_p_gdf)))
"""

# ╔═╡ f2597c01-7aa0-435e-a808-aee3c8a0e775
# ╠═╡ disabled = true
#=╠═╡
therm_p_gdf_iter[plot_iter]
  ╠═╡ =#

# ╔═╡ ea9f8f29-cebc-4ab0-82db-85e1fb5fe055
# ╠═╡ disabled = true
#=╠═╡
therm_e_gdf_iter[plot_iter]
  ╠═╡ =#

# ╔═╡ 1dd37d49-25aa-4eb8-bc96-52968ae302aa
md"""
Choose which frames to plot:
- Plasma frame: $plot_pf_binder
- Shock frame:  $plot_sf_binder
- ISM frame:    $plot_ISM_binder
"""

# ╔═╡ 6a790ec5-430f-4ad3-9315-dbda9b2a591f
md"""
Select which iteration to plot:
`plot_iter` = $(index_binder) (min: $(minimum(idx_therm_p_gdf)), max: $(maximum(idx_therm_p_gdf)))
"""

# ╔═╡ e916e00e-9b31-478b-b4b5-abca4259385e
# The primary plot in this notebook, i.e., one spectra, is through
# AlgebraOfGraphics, which is built on top of Makie.
let fig = Figure(), df = therm_p_gdf_iter[plot_iter]
    # Create a plotting specification, which specifies our data source, how to
    # transform our data (`map_layer`), and how to present the data (`visual_layer`).
    spec = data(df) * map_layer * (visual_layer + visual(Scatter; markersize))
    title = "dN/dp of thermal particles (protons), iteration $plot_iter"
    plt = draw!(fig[1, 1], spec, axis = (; title, axis_properties...))
    legend!(fig[1, 1], plt; legend_properties...)
    fig
end

# ╔═╡ 08b8d65c-bd5f-4267-9330-52d68aa4f7b5
let fig = Figure(), df = therm_e_gdf_iter[plot_iter]
    # Create a plotting specification, which specifies our data source, how to
    # transform our data (`map_layer`), and how to present the data (`visual_layer`).
    spec = data(df) * map_layer * (visual_layer + visual(Scatter; markersize))
    title = "dN/dp of thermal particles (electrons), iteration $plot_iter"
    plt = draw!(fig[1, 1], spec, axis = (; title, axis_properties...))
    legend!(fig[1, 1], plt; legend_properties...)
    fig
end

# ╔═╡ Cell order:
# ╟─a1e20290-55e9-4355-bc45-74c5f825565f
# ╠═b0d4ca90-bb62-11f0-80dd-b560b36acdf6
# ╠═63875dd0-d278-44c4-a85a-78ca6d7e9e10
# ╠═a08cb144-1980-46e6-b835-65c93fd5a265
# ╠═d3e5bcd8-31e0-424f-8626-5eb4b679b818
# ╠═b1833ac8-965e-45d5-bffd-8e31fac39939
# ╠═24564257-76f5-4200-b9a8-30180e853233
# ╠═e50e35b0-885a-4f34-b716-a0167f62d8f7
# ╠═9e779980-f661-45cd-bdd1-20e073df9ede
# ╟─d6dbd935-f720-464c-a3ed-9a8f26b8d287
# ╠═4f092e96-b696-4a79-a109-95c2ecee448f
# ╟─e368d071-20c3-433d-b884-c0ecb3887b5a
# ╟─537c5e04-31ef-4c5e-a80e-b255ed166ba1
# ╠═9b3120c2-2b26-41e5-a067-19cadd0c5136
# ╠═7bc31948-1d91-461d-8fa0-3a58aebeace4
# ╠═b6a0acf6-628c-40c5-b5f9-c47ba3814ab2
# ╠═667ae8fa-015e-4c30-b8b0-f1bc9646a5ab
# ╟─4e419982-323d-4cdf-b026-4ebec2cd0c39
# ╠═fa5460cd-c3e5-440c-8c13-94f8374fa643
# ╠═6764b6b6-63cf-4a7b-8750-65002ee0a728
# ╟─c8442393-614f-4bfa-aeab-c8fad1b3e0b8
# ╟─34fb2a5f-6c38-43a7-aec5-0c05d163ab4e
# ╟─3dd36273-6ab6-4003-a7f3-4ca4a456ed21
# ╟─76507571-215c-457f-8c5f-affdca267b1f
# ╟─01f0eec1-3f3e-44b3-b9cd-246eb3b9e557
# ╠═c14fef4c-6dcd-410c-8d23-8beee9a38e1e
# ╠═10e30fc9-7a25-4672-86e9-724d490a0db4
# ╟─0ce3bbd3-150f-4fb5-a7f6-4264b11dfbe6
# ╠═d5fe4e4c-c089-43fb-a8cc-26053ef5775f
# ╠═c4d125ac-7412-4f52-838c-33e1bfa23386
# ╠═25858c85-8006-4572-ad81-e8b220fab0d8
# ╠═4e8d3b62-0c32-427b-b1f2-e1a74d82a66c
# ╠═3b1ecd57-8f9f-40f7-818c-651a8a320f6f
# ╠═da63beaa-882c-4d34-b69a-3423d98e12c1
# ╟─78db62d2-4717-4bcf-a6fb-6cd80a2b9faf
# ╠═2e5e7272-0313-471d-8987-de64f79b84b0
# ╠═8c12b2c3-b5e3-4145-97b8-bba2ab5f03fc
# ╠═dc4d6726-25e1-4aeb-9029-fa88282667db
# ╠═6aae8a53-9e2a-4c54-b5da-af9452fa2d5b
# ╠═911db519-1f6c-4534-9de7-5e03e0869ae3
# ╟─720c2380-c92a-4a1f-a8f4-5e455a79651a
# ╟─6c256957-03b3-40e6-9961-e26568c2c72c
# ╟─f3e07ca1-3e4c-41fc-b877-dc635097b15e
# ╠═f2597c01-7aa0-435e-a808-aee3c8a0e775
# ╠═ea9f8f29-cebc-4ab0-82db-85e1fb5fe055
# ╟─1dd37d49-25aa-4eb8-bc96-52968ae302aa
# ╟─6a790ec5-430f-4ad3-9315-dbda9b2a591f
# ╠═e916e00e-9b31-478b-b4b5-abca4259385e
# ╠═08b8d65c-bd5f-4267-9330-52d68aa4f7b5
