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

# ╔═╡ 7a050dc5-7772-4933-959f-bf4fb478fc7d
using PlutoUI

# ╔═╡ 40efcd80-db38-4db3-a193-6e65ee5c4367
using PlutoUI: Slider

# ╔═╡ acbe8855-c586-46e3-a72b-a556df77b547
#using WGLMakie
using CairoMakie

# ╔═╡ 3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
using Missings

# ╔═╡ fe2b3846-c753-4685-8704-e6fb50624989
using Printf

# ╔═╡ a5526239-2f05-4618-8868-0f552855d574
md"""
# Preamble
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
    max-width: 80%;
    padding-left: max(360px, 10%);
    padding-right: 0%;
}
</style>"""

# ╔═╡ 8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
md"""
# Read data file
"""

# ╔═╡ ee3eab6d-7913-4650-a5c3-aabf0747a58a
const datadir = "G:/My Drive/MC Scattering/Processed-data";

# ╔═╡ 3bc899e5-ce26-4384-bf7b-f0bb3820f08d
CR_p_gdf_iter = load_object(joinpath(datadir, "dNdp-CR-protons-iteration-split.jld2"));

# ╔═╡ 80b16c96-b0f3-42a8-8544-7fbd9c06a1d9
CR_e_gdf_iter = load_object(joinpath(datadir, "dNdp-CR-electrons-iteration-split.jld2"));

# ╔═╡ 3a2e4aee-bc90-493e-84b1-79897934f16a
CR_p_gdf_iter

# ╔═╡ fae99b13-0b14-45c0-989f-8d0f22f0e96c
CR_e_gdf_iter

# ╔═╡ 628130bf-da25-4799-8e5e-3d2db15b1e49
md"""
# Plot Cosmic Ray data
"""

# ╔═╡ f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
const axis_properties = (xminorgridvisible = true, yminorgridvisible = true);

# ╔═╡ 50b1a87f-49ff-4d93-aa6e-f042a87b875e
const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ 3cc54622-c4e4-4c59-8828-4aa899a51e51
const markersize = 5;

# ╔═╡ 19a41e11-d031-498c-adbb-082e682fb67e
md"""
## Individual iterations
"""

# ╔═╡ ecf80697-b786-4b02-9563-f3d082383b76
md"""
Choose which frames to plot:
- Plasma frame: $(@bind do_plot_pf CheckBox(default=true))
- Shock frame: $(@bind do_plot_sf CheckBox(default=false))
- ISM frame: $(@bind do_plot_ISM CheckBox(default=false))
"""

# ╔═╡ 9c7fc2ae-f6b5-4a4d-90ed-63967ea55200
#plot_iter = 1

# ╔═╡ c2b3d96a-216e-4abe-8b0f-625419ac072f
#CR_e_gdf_iter[plot_iter]

# ╔═╡ 47a47a1d-4247-4a1a-a629-0a580253b41d
md"""
at plot_iter 1, there's a weird kick at the end. why?

Dr. Warren suggestion: momentum splitting. to be investigated
"""

# ╔═╡ d7d554cf-2f16-49e1-849d-25b5088e85ff
md"""
`plot_iter` = $(
    @bind plot_iter NumberField(axes(CR_p_gdf_iter, 1))
)
"""

# ╔═╡ 220c3ca5-e0b5-4f5c-86b0-e5d7cdd67558
let f = Figure(), df = CR_p_gdf_iter[plot_iter]
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (protons), iteration $plot_iter",
        xlabel = "log(p) (cgs)", ylabel = "log(dN/dp)", axis_properties...)

    p, N = df.log_p_nat, df.log_dNdp_cr_sf
    do_plot_sf && lines!(ax, p, N, label = "shock frame", color = color_sf_p)
    p, N = df.log_p_nat, df.log_dNdp_cr_pf
    do_plot_pf && lines!(ax, p, N, label = "plasma frame", color = color_pf_p)
    p, N = df.log_p_nat, df.log_dNdp_cr_ISM
    do_plot_ISM && lines!(ax, p, N, label = "ISM frame", color = color_ISM_p)

    #xlims!(ax, -16, -3)
    #ylims!(ax, 30, 60)
    axislegend(ax)
    f
end

# ╔═╡ 7879a41a-a284-452b-9505-a239209f1ed0
let f = Figure(), df = CR_e_gdf_iter[plot_iter]
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (electrons), iteration $plot_iter",
        xlabel = "log(p) (nat)", ylabel = "log(dN/dp)",
        axis_properties...)

    do_plot_sf && lines!(ax, df.log_p_nat, df.log_dNdp_cr_sf, label = "shock frame", color = color_sf_e)
    do_plot_pf && lines!(ax, df.log_p_nat, df.log_dNdp_cr_pf, label = "plasma frame", color = color_pf_e)
    do_plot_ISM && lines!(ax, df.log_p_nat, df.log_dNdp_cr_ISM, label = "ISM frame", color = color_ISM_e)

    #xlims!(ax, extrema(df.psd_mom_bounds_nat))
    axislegend(ax)
    f
end

# ╔═╡ 6537effb-12e6-4f4e-b34f-15dd33547921
const σ = 2.23;

# ╔═╡ f4930314-a64c-4b6a-bcef-c0d9dcf2ef81
let
    dfp = CR_p_gdf_iter[plot_iter]
    dfe = CR_e_gdf_iter[plot_iter]

    f = Figure()

    ax = Axis(
        f[1,1]; title = "dN/dp of Cosmic rays, iteration $plot_iter",
        xlabel = "log(p) (nat)", ylabel = "log dN/dp + σ log p", axis_properties...)

    if do_plot_pf
        p, N = (dfp.log_p_nat, dfp.log_dNdp_cr_pf)
        scatterlines!(ax, p, N .+ σ*p, label = "protons, plasma frame"; color = color_pf_p, markersize)

        p, N = (dfe.log_p_nat, dfe.log_dNdp_cr_pf)
        scatterlines!(ax, p, N .+ σ*p, label = "electrons, plasma frame"; color = color_pf_e, markersize)
    end
    if do_plot_sf
        p, N = (dfp.log_p_nat, dfp.log_dNdp_cr_sf)
        scatterlines!(ax, p, N .+ σ*p, label = "protons, shock frame"; color = color_sf_p, markersize)

        p, N = (dfe.log_p_nat, dfe.log_dNdp_cr_sf)
        scatterlines!(ax, p, N .+ σ*p, label = "electrons, shock frame"; color = color_sf_e, markersize)
    end
    if do_plot_ISM
        p, N = (dfp.log_p_nat, dfp.log_dNdp_cr_ISM)
        scatterlines!(ax, p, N .+ σ*p, label = "protons, ISM frame"; color = color_ISM_p, markersize)

        p, N = (dfe.log_p_nat, dfe.log_dNdp_cr_ISM)
        scatterlines!(ax, p, N .+ σ*p, label = "electrons, ISM frame"; color = color_ISM_e, markersize)
    end

    hlines!(ax, 57.8, color = color_pf_p)
    hlines!(ax, 56.5, color = color_pf_e)

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
## Multiple iterations
"""

# ╔═╡ 1f35f220-7739-4097-b51d-0ab6000be247
let f = Figure()
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (electrons)",
        xlabel = "log(p) (nat)", ylabel = "log(dN/dp) + σ log(p)",
        axis_properties...)

    for (i, dfe) in enumerate(CR_e_gdf_iter[5775:5779])
        p, N = dfe.log_p_nat, dfe.log_dNdp_cr_pf
        scatterlines!(ax, p, N + σ*p, label = "plasma frame (iter $i)"; markersize)
    end

    xlims!(ax, -0.3, 5)
    ylims!(ax, 56, 57.0)
    #axislegend(ax)
    f
end

# ╔═╡ 6c07e039-2575-49a6-a50d-531c40ee7965
let f = Figure()
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (protons)",
        xlabel = "log(p) (nat)", ylabel = "log dN/dp + σ log p",
        axis_properties...)

    for (i, dfproton) in enumerate(CR_p_gdf_iter[5620:5630])
        p, N = dfproton.log_p_nat, dfproton.log_dNdp_cr_pf
        scatterlines!(ax, p, N + σ*p, label = "plasma frame (iter $i)"; markersize)
    end

    #hlines!(ax, 57.5)

    xlims!(ax, 2, 5)
    ylims!(ax, 57.2, 58.3)
    #axislegend(ax, position = :lb)
    f
end

# ╔═╡ e5dbf380-3480-4d96-881a-8c562b5fc6ab
md"""
## All at once
"""

# ╔═╡ 4e26e9ec-b4f2-46f8-bada-945c00cb4907
# ╠═╡ disabled = true
#=╠═╡
let f = Figure()
    ax = Axis(
        f[1,1];
        title = "dN/dp of Cosmic rays (protons)",
        xlabel = "log(p) (nat)", ylabel = "log dN/dp + σ log p",
        axis_properties...)

    for (i, dfproton) in enumerate(CR_p_gdf_iter)
        p, N = dfproton.log_p_nat, dfproton.log_dNdp_cr_pf
        scatterlines!(ax, p, N + σ*p, label = "plasma frame (iter $i)"; markersize)
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
        p, N = dfe.log_p_nat, dfe.log_dNdp_cr_pf
        scatterlines!(ax, p, N + σ*p, label = "plasma frame (iter $i)"; markersize)
    end

    xlims!(ax, -0.3, 5)
    ylims!(ax, 56, 57.4)
    #axislegend(ax, position = :lb)
    f
end
  ╠═╡ =#

# ╔═╡ 8d03de5e-d344-4efd-b9af-dd5391028780
md"""
# Constants and functions
"""

# ╔═╡ 4272e1f1-1c3b-439f-9406-5315a901587b


# ╔═╡ Cell order:
# ╠═f1ee2cb0-8274-11ef-0826-f55183647219
# ╟─a5526239-2f05-4618-8868-0f552855d574
# ╟─cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
# ╠═7899ae97-fbc2-43e5-ac77-c6d725f0371e
# ╠═b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
# ╠═7a050dc5-7772-4933-959f-bf4fb478fc7d
# ╠═40efcd80-db38-4db3-a193-6e65ee5c4367
# ╠═acbe8855-c586-46e3-a72b-a556df77b547
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
# ╠═f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
# ╠═50b1a87f-49ff-4d93-aa6e-f042a87b875e
# ╠═3cc54622-c4e4-4c59-8828-4aa899a51e51
# ╟─19a41e11-d031-498c-adbb-082e682fb67e
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
# ╠═9c7fc2ae-f6b5-4a4d-90ed-63967ea55200
# ╠═c2b3d96a-216e-4abe-8b0f-625419ac072f
# ╟─47a47a1d-4247-4a1a-a629-0a580253b41d
# ╟─220c3ca5-e0b5-4f5c-86b0-e5d7cdd67558
# ╟─7879a41a-a284-452b-9505-a239209f1ed0
# ╟─d7d554cf-2f16-49e1-849d-25b5088e85ff
# ╟─f4930314-a64c-4b6a-bcef-c0d9dcf2ef81
# ╠═6537effb-12e6-4f4e-b34f-15dd33547921
# ╟─67f27108-eb9d-49b0-95ae-e016973e02b5
# ╟─1f35f220-7739-4097-b51d-0ab6000be247
# ╟─6c07e039-2575-49a6-a50d-531c40ee7965
# ╟─e5dbf380-3480-4d96-881a-8c562b5fc6ab
# ╠═4e26e9ec-b4f2-46f8-bada-945c00cb4907
# ╟─a0be5567-9256-4c03-9a96-11d4d1973347
# ╠═526cd197-9ec7-445a-9018-3163d3916e10
# ╟─8d03de5e-d344-4efd-b9af-dd5391028780
# ╠═4272e1f1-1c3b-439f-9406-5315a901587b
