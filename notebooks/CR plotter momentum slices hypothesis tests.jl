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
import Pkg; Pkg.activate(Base.current_project())

# ╔═╡ fb9581c0-15f2-4e26-9a1b-2f6ba2bcc5b9
using Revise

# ╔═╡ e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
using DrWatson: datadir

# ╔═╡ b0f9b523-5272-4735-8c0a-1694ffbf72fd
using JLD2: load_object

# ╔═╡ 7899ae97-fbc2-43e5-ac77-c6d725f0371e
using DataFrames

# ╔═╡ b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
using Distributions

# ╔═╡ 800673a1-dcb4-471a-a628-74a92aee2941
using MCScatteringDataAnalysis

# ╔═╡ d8a66ccd-efdc-4b85-9af3-cfa5624c88e8
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

# ╔═╡ 49902e99-870d-4d19-afb0-1de612c185df
using StatsBase

# ╔═╡ 4ac32bef-af81-4f7e-8e97-7eac4dd2bf69
using LinearAlgebra

# ╔═╡ fd47dab7-426c-44fc-8038-00e378324e41
using HypothesisTests

# ╔═╡ f0e77bbd-e420-49f1-9b40-f9d994888b93
md"""
# Plot flux hypothesis tests for each momentum slice
"""

# ╔═╡ cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
md"""
## Import packages
"""

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

# ╔═╡ 802ee39c-3d64-40eb-b7e7-0c032ff52c84
CR_p_gdf_momentum_filename = datadir("Lorentz-5-processed", "dNdp-CR-protons-momentum-split.jld2");

# ╔═╡ 754ec337-6b5e-4606-ac1f-4ba8ad251631
CR_e_gdf_momentum_filename = datadir("Lorentz-5-processed", "dNdp-CR-electrons-momentum-split.jld2");

# ╔═╡ bdb9591b-b7ac-47e6-98bc-f18921bb64f9
CR_p_gdf_momentum = load_object(CR_p_gdf_momentum_filename);

# ╔═╡ 3777306e-eb41-413b-80a9-72cdc0228a94
CR_e_gdf_momentum = load_object(CR_e_gdf_momentum_filename);

# ╔═╡ 68c8329f-501e-47df-8047-d3cbc319e705
md"""
For protons: (enable the next cell to view Grouped-DataFrame statistics)
"""

# ╔═╡ a36ea9cf-176f-40bd-8577-cc2ea8db64af
# ╠═╡ disabled = true
#=╠═╡
CR_gdfstats(CR_p_gdf_momentum)
  ╠═╡ =#

# ╔═╡ 985a2460-3fbc-4935-af59-2e734786c973
md"""
For electrons: (enable the next cell to view Grouped-DataFrame statistics)
"""

# ╔═╡ d85427f4-86ed-4c04-980a-a4152b5875e8
# ╠═╡ disabled = true
#=╠═╡
CR_gdfstats(CR_e_gdf_momentum)
  ╠═╡ =#

# ╔═╡ 1572a05b-77db-43a4-81cd-d9eb3c9bf2e0
idx_CR_p_gdf = axes(CR_p_gdf_momentum, 1);

# ╔═╡ f8d26a2c-2789-4b04-8bb7-71bf19686bbd
idx_CR_e_gdf = axes(CR_e_gdf_momentum, 1);

# ╔═╡ 5c65e63d-e6bd-41fa-adbf-0e1717075956
md"""
UI element variables for selecting which momentum slice we want to inspect
"""

# ╔═╡ b679ff7a-89f4-4f92-9dec-2ba3c538d715
proton_index_binder = @bind proton_momentum_index NumberField(idx_CR_p_gdf, default = 13);

# ╔═╡ 59f400ef-3b5b-424c-ae08-41b94d220ce9
electron_index_binder = @bind electron_momentum_index NumberField(idx_CR_e_gdf, default = 13);

# ╔═╡ 628130bf-da25-4799-8e5e-3d2db15b1e49
md"""
# Plot Cosmic Ray data
"""

# ╔═╡ ecf233ad-d75e-4aa5-bf7e-ff3e7b1d8755
md"""
## Plotting configuration
"""

# ╔═╡ 59a22149-3397-4e97-9f7b-5d502aacf293
markersize = 5

# ╔═╡ f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
axis_properties = (
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorticksvisible = true,
    yminorticksvisible = true,
    xlabel = "log(dN/dp)",
)

# ╔═╡ f86707a1-9d79-4df8-8798-3f7ea1d1797c


# ╔═╡ 50b1a87f-49ff-4d93-aa6e-f042a87b875e
color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ 1b9f507c-1585-4ad1-8090-bdde6de972d6
md"""
Define the PlutoUI binders for selecting which frames to plot
"""

# ╔═╡ d6513272-7232-43e6-ac88-a58462181041
plot_pf_binder = @bind do_plot_pf CheckBox(default = true);

# ╔═╡ 59be6983-6e37-4a70-8929-69176a5f807e
plot_sf_binder = @bind do_plot_sf CheckBox(default = false);

# ╔═╡ 60deb76f-3efe-4e0d-b176-9f0169259dca
plot_ISM_binder = @bind do_plot_ISM CheckBox(default = false);

# ╔═╡ 105361e9-cafd-4755-bcbd-fdcbcb07b291
map_layer = let
    x_map = :log_dNdp_cr_pf => "log(dN/dp)"
    y_label = "log(dN/dp)"

    pf_map = mapping(x_map, :log_dNdp_cr_pf => y_label, color = direct("plasma frame"))
    sf_map = mapping(x_map, :log_dNdp_cr_sf => y_label, color = direct("shock frame"))
    ISM_map = mapping(x_map, :log_dNdp_cr_ISM => y_label, color = direct("ISM frame"))

    combined_layer = zerolayer()
    if do_plot_pf
        combined_layer += pf_map
    end
    if do_plot_sf
        combined_layer += sf_map
    end
    if do_plot_ISM
        combined_layer += ISM_map
    end
    combined_layer
end

# ╔═╡ 718f3f64-b581-41b2-86aa-d561cc71e6f8
plot_electrons_binder = @bind do_plot_electrons CheckBox(default = false);

# ╔═╡ ce8b1307-dc78-463b-9f41-04fe5dded525
md"""
## Histograms
"""

# ╔═╡ ecf80697-b786-4b02-9563-f3d082383b76
md"""
Choose which frames to plot:
- Plasma frame: $plot_pf_binder
- Shock frame: $plot_sf_binder
- ISM frame: $plot_ISM_binder
"""

# ╔═╡ f3132403-113d-4b30-9fd0-379d28ade3c7
md"""
# Normal distribution inference
"""

# ╔═╡ 8bc20375-2562-4611-a67b-5884aa99b5f0
md"""
## Statistical estimates (MLE/MoM/...)
"""

# ╔═╡ e6b9701d-3d27-4c0c-b0b9-9879527f369c
normal_distrib_protons = fitdistributions(v -> fitdistribution(Normal, v), CR_p_gdf_momentum)
# normal_distrib_protons = fitdistributions(fitnormal, CR_p_gdf_momentum)

# ╔═╡ e75ea9c0-59ca-4097-b4f6-6a3af04dc308
normal_distrib_electrons = fitdistributions(v -> fitdistribution(Normal, v), CR_e_gdf_momentum)

# ╔═╡ da107273-c428-4c68-80a9-8f82cb211497
md"""
# Hypothesis tests
"""

# ╔═╡ f330af91-60a6-46ac-bdc5-ec49c216fccb
md"""
Get goodness of fits
"""

# ╔═╡ d77a95bf-2d54-46d5-81ca-b671ee1db695
p_values_scale_checkbox_binder = @bind plot_p_values_in_logscale CheckBox();

# ╔═╡ 04dad413-0dc0-4ceb-81c2-e208ef082f38
p_val_yscale = plot_p_values_in_logscale ? log10 : identity;

# ╔═╡ 98675d19-3b1b-4be0-9e48-ab0ffd019647
md"""
## Anderson–Darling test
"""

# ╔═╡ 2e79471f-3430-4b1c-91fe-80434de63cb2
ad_scores_p = get_ad_scores(CR_p_gdf_momentum, normal_distrib_protons.pf, col = :log_dNdp_cr_pf);

# ╔═╡ bd8f636c-6033-434e-a220-a07397679431
ad_scores_e = get_ad_scores(CR_e_gdf_momentum, normal_distrib_electrons.pf, col = :log_dNdp_cr_pf);

# ╔═╡ 89b7b4d0-fb46-4436-a48d-0731cef1dc2c
md"""
Plot p-values in log scale? (Uncheck for linear)
$p_values_scale_checkbox_binder
"""

# ╔═╡ e68032c7-36bd-4817-9bae-09e7dcb04802
md"""
Should we plot electrons? $plot_electrons_binder
"""

# ╔═╡ b499bf86-3e7a-441a-809a-934a1a8dd402
md"""
## Shapiro–Wilk test
"""

# ╔═╡ a2dca585-2b84-4958-8ba6-af51602c4d8a
sw_scores_p = get_sw_scores(CR_p_gdf_momentum, col = :log_dNdp_cr_pf);

# ╔═╡ ec6883c9-b6bb-4e7e-bd6d-e65d6e06144d
sw_scores_e = get_sw_scores(CR_e_gdf_momentum, col = :log_dNdp_cr_pf);

# ╔═╡ 1b8925de-a9fe-4326-ad1a-a0a022ccbc6b
md"""
Plot p-values in log scale? (Uncheck for linear)
$p_values_scale_checkbox_binder
"""

# ╔═╡ 79736549-5c2c-4193-9757-16a57b0a535d
md"""
Should we plot electrons? $plot_electrons_binder
"""

# ╔═╡ 2ab2979f-1ad4-4168-b59c-a25e57d4826a
md"""
## Kolmogorov–Smirnov test
"""

# ╔═╡ fc4dddd0-cfca-407e-ad94-622f53b148b3
ks_scores_p = get_ks_scores(CR_p_gdf_momentum, normal_distrib_protons.pf, col = :log_dNdp_cr_pf);

# ╔═╡ cfdbbdea-c4fe-44ac-9de5-70720a138286
ks_scores_e = get_ks_scores(CR_e_gdf_momentum, normal_distrib_electrons.pf, col = :log_dNdp_cr_pf);

# ╔═╡ b825855d-bdfc-4aaa-ad9e-82ee4e8d1201
md"""
Plot p-values in log scale? (Uncheck for linear)
$p_values_scale_checkbox_binder
"""

# ╔═╡ cabc5f90-2e66-41fc-956a-8d2cd0b36bf5
md"""
Should we plot electrons? $plot_electrons_binder
"""

# ╔═╡ 8d03de5e-d344-4efd-b9af-dd5391028780
md"""
# Constants and functions
"""

# ╔═╡ 377aaf8f-b909-4c42-bc77-912fd300c300
normalization = :pdf;

# ╔═╡ 32edc221-e586-4510-9427-977b22f62f6c
md"""
Vector of momentum slices
"""

# ╔═╡ e8406a6a-ecc2-49d2-b67a-503b4ef5764b
proton_log_p_nat = keys(CR_p_gdf_momentum) .|> values .|> first;

# ╔═╡ 71404de8-f8b2-4d26-b7d7-41064cae1447
log_p_nat_at_slice_p = proton_log_p_nat[proton_momentum_index];

# ╔═╡ 89bcb29b-0b1c-4e3a-91cb-282c05df2bc5
md"""
Value of proton momentum at current slice: log(*p*/*m*ₚ*c*) = $log_p_nat_at_slice_p
"""

# ╔═╡ 589661b1-6a64-4db5-ac40-c1565c29c3cc
electron_log_p_nat = keys(CR_e_gdf_momentum) .|> values .|> first;

# ╔═╡ 6c16fc5a-7113-4b6e-abf2-de1275cceda5
log_p_nat_at_slice_e = electron_log_p_nat[electron_momentum_index];

# ╔═╡ c9b9969c-2c7f-436e-b5a1-603138a4e196
md"""
Value of electron momentum at current slice: log(*p*/*m*ₚ*c*) = $log_p_nat_at_slice_e
"""

# ╔═╡ cee91c99-adc0-4185-a7c3-e2164b95a003
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Anderson–Darling p-value vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)",
        yscale = p_val_yscale,
    )

    scatterlines!(ax, proton_log_p_nat, passmissing(pvalue).(ad_scores_p), color = color_pf_p, label = "protons, plasma frame"; markersize)
    do_plot_electrons && scatterlines!(ax, electron_log_p_nat, passmissing(pvalue).(ad_scores_e), color = color_pf_e, label = "electrons, plasma frame"; markersize)

    axislegend(ax, position = :lt, framevisible = false)

    fig
end

# ╔═╡ a49ff5ab-6077-4bb2-b694-6f3662982745
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Shapiro–Wilk p-value vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)",
        yscale = p_val_yscale,
    )

    scatterlines!(ax, proton_log_p_nat, passmissing(pvalue).(sw_scores_p), color = color_pf_p, label = "protons, plasma frame"; markersize)
    do_plot_electrons && scatterlines!(ax, electron_log_p_nat, passmissing(pvalue).(sw_scores_e), color = color_pf_e, label = "electrons, plasma frame"; markersize)

    axislegend(ax, position = plot_p_values_in_logscale ? :cb : :lt)

    fig
end

# ╔═╡ 08542eea-964a-4f1d-aae5-2b50a628588a
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Kolmogorov–Smirnov p-value vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)",
        yscale = p_val_yscale,
    )

    scatterlines!(ax, proton_log_p_nat, passmissing(pvalue).(ks_scores_p), color = color_pf_p, label = "protons, plasma frame"; markersize)
    do_plot_electrons && scatterlines!(ax, electron_log_p_nat, passmissing(pvalue).(ks_scores_e), color = color_pf_e, label = "electrons, plasma frame"; markersize)

    axislegend(ax, position = plot_p_values_in_logscale ? :cb : :lt, framevisible = false)

    fig
end

# ╔═╡ Cell order:
# ╟─f0e77bbd-e420-49f1-9b40-f9d994888b93
# ╟─cd809ca8-2cc4-435d-ab8b-b7b24fa40ed1
# ╠═f1ee2cb0-8274-11ef-0826-f55183647219
# ╠═e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
# ╠═fb9581c0-15f2-4e26-9a1b-2f6ba2bcc5b9
# ╠═b0f9b523-5272-4735-8c0a-1694ffbf72fd
# ╠═7899ae97-fbc2-43e5-ac77-c6d725f0371e
# ╠═b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
# ╠═800673a1-dcb4-471a-a628-74a92aee2941
# ╠═d8a66ccd-efdc-4b85-9af3-cfa5624c88e8
# ╠═547aad6f-32db-405d-9886-a727f1591101
# ╠═7a050dc5-7772-4933-959f-bf4fb478fc7d
# ╠═3791e767-dcf1-4f9d-909d-a7d08e4c5f9c
# ╠═fe2b3846-c753-4685-8704-e6fb50624989
# ╟─5c6b130f-0a51-4131-bab7-40b059c4cc11
# ╠═b544df91-fe2d-4396-892c-7faea2edd141
# ╟─4415022a-54dc-4f3d-a651-f66ae63dd051
# ╟─8dfe6f3c-f693-4c73-8152-8c43c1c1ff42
# ╠═802ee39c-3d64-40eb-b7e7-0c032ff52c84
# ╠═754ec337-6b5e-4606-ac1f-4ba8ad251631
# ╠═bdb9591b-b7ac-47e6-98bc-f18921bb64f9
# ╠═3777306e-eb41-413b-80a9-72cdc0228a94
# ╟─68c8329f-501e-47df-8047-d3cbc319e705
# ╠═a36ea9cf-176f-40bd-8577-cc2ea8db64af
# ╟─985a2460-3fbc-4935-af59-2e734786c973
# ╠═d85427f4-86ed-4c04-980a-a4152b5875e8
# ╠═1572a05b-77db-43a4-81cd-d9eb3c9bf2e0
# ╠═f8d26a2c-2789-4b04-8bb7-71bf19686bbd
# ╟─5c65e63d-e6bd-41fa-adbf-0e1717075956
# ╠═b679ff7a-89f4-4f92-9dec-2ba3c538d715
# ╠═59f400ef-3b5b-424c-ae08-41b94d220ce9
# ╟─628130bf-da25-4799-8e5e-3d2db15b1e49
# ╟─ecf233ad-d75e-4aa5-bf7e-ff3e7b1d8755
# ╟─59a22149-3397-4e97-9f7b-5d502aacf293
# ╟─f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
# ╠═f86707a1-9d79-4df8-8798-3f7ea1d1797c
# ╠═50b1a87f-49ff-4d93-aa6e-f042a87b875e
# ╟─105361e9-cafd-4755-bcbd-fdcbcb07b291
# ╟─1b9f507c-1585-4ad1-8090-bdde6de972d6
# ╠═d6513272-7232-43e6-ac88-a58462181041
# ╠═59be6983-6e37-4a70-8929-69176a5f807e
# ╠═60deb76f-3efe-4e0d-b176-9f0169259dca
# ╠═718f3f64-b581-41b2-86aa-d561cc71e6f8
# ╟─ce8b1307-dc78-463b-9f41-04fe5dded525
# ╠═71404de8-f8b2-4d26-b7d7-41064cae1447
# ╠═6c16fc5a-7113-4b6e-abf2-de1275cceda5
# ╟─89bcb29b-0b1c-4e3a-91cb-282c05df2bc5
# ╟─c9b9969c-2c7f-436e-b5a1-603138a4e196
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
# ╟─f3132403-113d-4b30-9fd0-379d28ade3c7
# ╟─8bc20375-2562-4611-a67b-5884aa99b5f0
# ╠═e6b9701d-3d27-4c0c-b0b9-9879527f369c
# ╠═e75ea9c0-59ca-4097-b4f6-6a3af04dc308
# ╟─da107273-c428-4c68-80a9-8f82cb211497
# ╠═49902e99-870d-4d19-afb0-1de612c185df
# ╠═4ac32bef-af81-4f7e-8e97-7eac4dd2bf69
# ╟─f330af91-60a6-46ac-bdc5-ec49c216fccb
# ╠═fd47dab7-426c-44fc-8038-00e378324e41
# ╠═d77a95bf-2d54-46d5-81ca-b671ee1db695
# ╟─04dad413-0dc0-4ceb-81c2-e208ef082f38
# ╟─98675d19-3b1b-4be0-9e48-ab0ffd019647
# ╠═2e79471f-3430-4b1c-91fe-80434de63cb2
# ╠═bd8f636c-6033-434e-a220-a07397679431
# ╟─89b7b4d0-fb46-4436-a48d-0731cef1dc2c
# ╟─e68032c7-36bd-4817-9bae-09e7dcb04802
# ╟─cee91c99-adc0-4185-a7c3-e2164b95a003
# ╟─b499bf86-3e7a-441a-809a-934a1a8dd402
# ╠═a2dca585-2b84-4958-8ba6-af51602c4d8a
# ╠═ec6883c9-b6bb-4e7e-bd6d-e65d6e06144d
# ╟─1b8925de-a9fe-4326-ad1a-a0a022ccbc6b
# ╟─79736549-5c2c-4193-9757-16a57b0a535d
# ╟─a49ff5ab-6077-4bb2-b694-6f3662982745
# ╟─2ab2979f-1ad4-4168-b59c-a25e57d4826a
# ╠═fc4dddd0-cfca-407e-ad94-622f53b148b3
# ╠═cfdbbdea-c4fe-44ac-9de5-70720a138286
# ╟─b825855d-bdfc-4aaa-ad9e-82ee4e8d1201
# ╟─cabc5f90-2e66-41fc-956a-8d2cd0b36bf5
# ╟─08542eea-964a-4f1d-aae5-2b50a628588a
# ╟─8d03de5e-d344-4efd-b9af-dd5391028780
# ╠═377aaf8f-b909-4c42-bc77-912fd300c300
# ╟─32edc221-e586-4510-9427-977b22f62f6c
# ╠═e8406a6a-ecc2-49d2-b67a-503b4ef5764b
# ╠═589661b1-6a64-4db5-ac40-c1565c29c3cc
