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

# ╔═╡ 7899ae97-fbc2-43e5-ac77-c6d725f0371e
using JLD2, DataFrames

# ╔═╡ b137e7fa-f2ce-4cb1-85d7-87078a9aa9cc
using Distributions

# ╔═╡ 800673a1-dcb4-471a-a628-74a92aee2941
using MCScatteringDataAnalysis

# ╔═╡ d8a66ccd-efdc-4b85-9af3-cfa5624c88e8
# using CairoMakie
using WGLMakie
# using GLMakie

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

# ╔═╡ 70f93b37-a977-4dce-9fdd-a0497603a864
using MCScatteringDataAnalysis: get_onesample_scores

# ╔═╡ f0e77bbd-e420-49f1-9b40-f9d994888b93
md"""
# Plot fluxes for each momentum slice
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
For protons:
"""

# ╔═╡ a36ea9cf-176f-40bd-8577-cc2ea8db64af
CR_gdfstats(CR_p_gdf_momentum)

# ╔═╡ 985a2460-3fbc-4935-af59-2e734786c973
md"""
For electrons:
"""

# ╔═╡ d85427f4-86ed-4c04-980a-a4152b5875e8
CR_gdfstats(CR_e_gdf_momentum)

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
bins = 90;

# ╔═╡ 50b1a87f-49ff-4d93-aa6e-f042a87b875e
color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();

# ╔═╡ 1b9f507c-1585-4ad1-8090-bdde6de972d6
md"""
Define the PlutoUI binders for selecting which frames to plot
"""

# ╔═╡ d6513272-7232-43e6-ac88-a58462181041
plot_pf_binder = @bind do_plot_pf CheckBox(default=true);

# ╔═╡ 59be6983-6e37-4a70-8929-69176a5f807e
plot_sf_binder = @bind do_plot_sf CheckBox(default=false);

# ╔═╡ 60deb76f-3efe-4e0d-b176-9f0169259dca
plot_ISM_binder = @bind do_plot_ISM CheckBox(default=false);

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

# ╔═╡ 4553a97d-6b78-4268-90de-d8bee348d3d4
plot_electrons_binder = @bind do_plot_electrons CheckBox(default=true);

# ╔═╡ 3bf64608-0fa2-4fcb-9782-fd7a8de47bda
md"""
## Sample statistics
"""

# ╔═╡ 70717f68-e97b-401e-bcf7-0684ade30b07
gdf_sample_stats(statistic, gdf; column) = [
    statistic(collect(skipmissing(df[!,column]))) for df in gdf];

# ╔═╡ c09cf00d-0a07-4159-9009-45afdb8343fb
CR_p_mean_log_dNdp = (;
    pf = gdf_sample_stats(mean, CR_p_gdf_momentum; column = :log_dNdp_cr_pf),
    sf = gdf_sample_stats(mean, CR_p_gdf_momentum; column = :log_dNdp_cr_sf),
    ISM = gdf_sample_stats(mean, CR_p_gdf_momentum; column = :log_dNdp_cr_ISM),
);

# ╔═╡ a8fdec46-2d8d-4a06-800a-7e03c64468dd
CR_e_mean_log_dNdp = (;
    pf = gdf_sample_stats(mean, CR_e_gdf_momentum; column = :log_dNdp_cr_pf),
    sf = gdf_sample_stats(mean, CR_e_gdf_momentum; column = :log_dNdp_cr_sf),
    ISM = gdf_sample_stats(mean, CR_e_gdf_momentum; column = :log_dNdp_cr_ISM),
);

# ╔═╡ a9b9b9aa-b850-49de-8412-2930e7004a36
CR_p_std_log_dNdp = (;
    pf = gdf_sample_stats(std, CR_p_gdf_momentum; column = :log_dNdp_cr_pf),
    sf = gdf_sample_stats(std, CR_p_gdf_momentum; column = :log_dNdp_cr_sf),
    ISM = gdf_sample_stats(std, CR_p_gdf_momentum; column = :log_dNdp_cr_ISM),
);

# ╔═╡ 734626db-7646-4a0c-9303-2aec6e371159
CR_e_std_log_dNdp = (;
    pf = gdf_sample_stats(std, CR_e_gdf_momentum; column = :log_dNdp_cr_pf),
    sf = gdf_sample_stats(std, CR_e_gdf_momentum; column = :log_dNdp_cr_sf),
    ISM = gdf_sample_stats(std, CR_e_gdf_momentum; column = :log_dNdp_cr_ISM),
);

# ╔═╡ e85bc3ef-97a4-4d0e-937e-3ca290a28086
md"""
Choose which frames to plot:
- Plasma frame: $plot_pf_binder
- Shock frame: $plot_sf_binder
- ISM frame: $plot_ISM_binder
"""

# ╔═╡ 88527827-9710-4c63-964f-691aa8909d1f
md"""
Should we plot electrons? $plot_electrons_binder
"""

# ╔═╡ cd684b62-1aad-404a-87b4-e7b406c0c989
md"""
All of the plotting code below can be simplified much, much more with AlgebraOfGraphics. The current code is a travesty. :(
"""

# ╔═╡ 54ec7213-3e48-423b-9208-befc6583b908
md"""
### Sample size
"""

# ╔═╡ 932c2a77-0198-4df4-a4bd-30d0bda93946
md"""
### Means
"""

# ╔═╡ 7da22399-718c-4f16-ba02-7ae27773942b
md"""
#### Mean with uncertainty envelope
"""

# ╔═╡ 5767b9ac-64c2-4d2f-ad42-961184c7edc7
md"""
### Standard deviations
"""

# ╔═╡ 7495e7e9-3d50-4401-baef-d2e3c11e6b46
md"""
### Skewness
"""

# ╔═╡ bb3a74ce-78ee-487e-a413-4c0e035e8818
md"""
### Kurtosis
"""

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

# ╔═╡ febdc8a1-00bb-47a7-83d2-6cccef5190f5
# ╠═╡ disabled = true
#=╠═╡
CR_p_gdf_momentum[proton_momentum_index]
  ╠═╡ =#

# ╔═╡ 35710ad9-f2e4-487b-be19-c29500633726
md"""
Proton momentum slice to plot (index): $proton_index_binder (min: $(minimum(idx_CR_p_gdf)), max: $(maximum(idx_CR_p_gdf)))
"""



# ╔═╡ 7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
md"""
Electron momentum slice to plot (index): $electron_index_binder (min: $(minimum(idx_CR_e_gdf)), max: $(maximum(idx_CR_e_gdf)))
"""

# ╔═╡ 9ea7a3a4-987d-416d-88d1-672e3cce23c5
md"""
## dN/dp vs. iteration
"""

# ╔═╡ b7a96870-784e-4ce0-830d-d245fc16e5f4
# ╠═╡ disabled = true
#=╠═╡
let df = CR_p_gdf_momentum[proton_momentum_index]

    fig = Figure()
    ax = Axis(
        fig[1,1];
        title = "dN/dp of Cosmic rays (protons) against iteration, momentum slice",
        xlabel = "Iteration", ylabel = "log(dN/dp)")

    do_plot_pf && scatter!(ax, df.iter, df.log_dNdp_cr_pf, label = "plasma frame"; markersize, color = color_pf_p)
    do_plot_sf && scatter!(ax, df.iter, df.log_dNdp_cr_sf, label = "shock frame"; markersize, color = color_sf_p)
    do_plot_ISM && scatter!(ax, df.iter, df.log_dNdp_cr_ISM, label = "ISM frame"; markersize, color = color_ISM_p)

    #xlims!(ax, -16, -3)
    #ylims!(ax, -100, -98)
    leg = axislegend(ax, position = :rb, framevisible = false)
    #Legend(fig[1,2], ax)
    fig
end
  ╠═╡ =#

# ╔═╡ 4ac1798d-ec27-4571-9b2a-44cb432ef0d6
# ╠═╡ disabled = true
#=╠═╡
let df = CR_p_gdf_momentum[electron_momentum_index]

    fig = Figure()
    ax = Axis(
        fig[1,1];
        #aspect = AxisAspect(1.2),
        title = "dN/dp of Cosmic rays (electrons) against iteration, momentum slice",
        #axis_properties...,
        xlabel = "Iteration", ylabel = "log(dN/dp)")

    do_plot_pf && scatter!(ax, df.iter, df.log_dNdp_cr_pf, label = "plasma frame"; markersize, color = color_pf_e)
    do_plot_sf && scatter!(ax, df.iter, df.log_dNdp_cr_sf, label = "shock frame"; markersize, color = color_sf_e)
    do_plot_ISM && scatter!(ax, df.iter, df.log_dNdp_cr_ISM, label = "ISM frame"; markersize, color = color_ISM_e)

    #xlims!(ax, -16, -3)
    #ylims!(ax, -100, -98)
    axislegend(ax, position = :rb, framevisible = false)
    #Legend(fig[1,2], ax)
    fig
end
  ╠═╡ =#

# ╔═╡ a6a63cb1-1a13-4cc2-9730-b78dd3d3aee4
md"""
# Inspect tail
"""

# ╔═╡ 5dc367ca-2882-4b98-8f29-2b5390426a9b
log_dNdp = CR_p_gdf_momentum[proton_momentum_index].log_dNdp_cr_pf |> skipmissing |> collect;

# ╔═╡ 464e92d7-e414-4ddd-a81e-978f271961b2
md"""
Proton momentum slice to plot (index): $proton_index_binder (min: $(minimum(idx_CR_p_gdf)), max: $(maximum(idx_CR_p_gdf)))
"""

# ╔═╡ daded5be-4319-4a0a-8491-84a7273a844b
md"""
Sum-of-squared errors for analytical fit and histogram curve fit:
"""

# ╔═╡ aebbf1a7-d047-4fa6-aa36-4b9ae8b68127
md"""
(log-)Likelihoods for analytical fit and histogram curve fit:
"""

# ╔═╡ 29ec59ad-0e22-462a-ab6d-2065a56fc001
x, y = get_hist_curve(log_dNdp; nbins=bins)

# ╔═╡ 6cb898b3-98c5-4f3a-8d77-3deef7cf5358
md"""
---
"""

# ╔═╡ 2374b968-1172-48db-8ddd-7b4deae7817c
md"""
Truncate all middle data
"""

# ╔═╡ c79d7166-e0a5-465e-b78d-6eeee331a99c
cutoff_log_dNdp = (; low = 31, high = 33.8);

# ╔═╡ a2c5161c-7cbd-4314-8fda-59a8e1750da1
log_dNdp_cur_lowtail = filter(x -> x ≤ cutoff_log_dNdp.low, log_dNdp);

# ╔═╡ e5d14c88-fea2-4117-97ab-0aebf3711a5b
log_dNdp_cur_hightail = filter(x -> x ≥ cutoff_log_dNdp.high, log_dNdp);

# ╔═╡ 59444b54-893e-4f4e-b746-97de78417043
log_dNdp_cur_trunc = filter(x -> cutoff_log_dNdp.low ≤ x ≤ cutoff_log_dNdp.high, log_dNdp);

# ╔═╡ 316b7ed7-4fe8-40dc-8b1a-551a0100c57c
log_dNdp_cur_lowtail |> length

# ╔═╡ 07300c09-2361-40e0-a502-b018496184c8
log_dNdp_cur_hightail |> length

# ╔═╡ afabc297-408f-4643-8296-40be885adafc
log_dNdp_cur_trunc |> length

# ╔═╡ 60bd9873-0246-432d-9e68-bfe2aa0956b2
log_dNdp |> length

# ╔═╡ 46607c86-7310-4d81-9816-2283d28d1420
md"""
## Plot of tail
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

# ╔═╡ 54452e38-227e-4d06-ae74-7347aae2c021
fitted_dist_MLE = normal_distrib_protons.pf[proton_momentum_index]

# ╔═╡ 5bbd6e99-87e1-401c-a09e-065e2d426370
SSE_hist(log_dNdp, fitted_dist_MLE)

# ╔═╡ 55d8c831-27e6-4914-a836-7a05281e8fb3
sum(logpdf.(fitted_dist_MLE, log_dNdp))

# ╔═╡ e75ea9c0-59ca-4097-b4f6-6a3af04dc308
normal_distrib_electrons = fitdistributions(v -> fitdistribution(Normal, v), CR_e_gdf_momentum)

# ╔═╡ 639ab710-c411-4831-9e1d-d7fba723b7bc
md"""
## Curve-fitting to histogram
"""

# ╔═╡ 3860d0cf-20f4-4256-9286-8757afc38ef9
md"""
Approach it like curve-fitting
"""

# ╔═╡ f2fe3844-2be8-4da6-9656-40312304556b
normal_distrib_protons_from_curves = fitdistributions(v -> fit_dist_to_histogram(Normal, v; nbins=bins), CR_p_gdf_momentum)

# ╔═╡ b0d555b3-5087-4405-8343-ce304d482ca9
fitted_dist_curve = normal_distrib_protons_from_curves.pf[proton_momentum_index]

# ╔═╡ 5ab05dc9-3a98-4297-a47b-c4e0111b8c51
SSE_hist(log_dNdp, fitted_dist_curve)

# ╔═╡ 89f8d7a8-ea2e-4906-9460-da16154b0404
sum(logpdf.(fitted_dist_curve, log_dNdp))

# ╔═╡ 97291776-74f0-428a-ab4f-3c498b630000
normal_distrib_electrons_from_curves = fitdistributions(v -> fit_dist_to_histogram(Normal, v; nbins=bins), CR_e_gdf_momentum)

# ╔═╡ 452c9b2f-7138-4310-b0c6-df2be7ab8c76
md"""
## Distribution agreement curve
"""

# ╔═╡ b238afe1-3d1f-4e15-bc49-1b015a39c02c
proton_distances = bcdistances(normal_distrib_protons.pf, normal_distrib_protons_from_curves.pf)

# ╔═╡ 78a22648-c76a-4b5c-b552-9be000a60109
electron_distances = bcdistances(normal_distrib_electrons.pf, normal_distrib_electrons_from_curves.pf)

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

# ╔═╡ 94a91acd-a878-4c3c-9716-8bed60bf8c6c
md"""
## Root-sum-squared errors
"""

# ╔═╡ 222df0cb-0760-48a2-902e-91d32e451a11
sse_scores_p = get_sse_scores(CR_p_gdf_momentum, normal_distrib_protons.pf, col = :log_dNdp_cr_pf)

# ╔═╡ 32f07cd2-f62f-41e0-9211-8ac333bdd98d
sse_scores_p |> skipmissing |> findmax

# ╔═╡ c5d56d28-15d5-464f-9055-7bbfc9826e72
# sse_scores_p_curve = get_sse_scores(CR_p_gdf_momentum, normal_distrib_protons_from_curves.pf, col = :log_dNdp_cr_pf)
sse_scores_p_curve = get_onesample_scores((t...; kwd...) -> SSE_hist(t...; kwd..., relative = true, bias = 10), CR_p_gdf_momentum, normal_distrib_protons_from_curves.pf; col = :log_dNdp_cr_pf)

# ╔═╡ cbea4ff4-b132-4abb-97c6-e406a339ced6
sse_scores_e = get_sse_scores(CR_e_gdf_momentum, normal_distrib_electrons.pf, col = :log_dNdp_cr_pf)

# ╔═╡ 85feafa4-a572-40ca-9975-fb0d3d5309f7
sse_scores_e_curve = get_onesample_scores((t...; kwd...) -> SSE_hist(t...; kwd..., relative = true, bias = 0.5), CR_e_gdf_momentum, normal_distrib_electrons_from_curves.pf; col = :log_dNdp_cr_pf)

# ╔═╡ 79dc57bb-d66d-4608-a775-9dfc58af1995
md"""
Plot p-values in log scale? (Uncheck for linear)
$p_values_scale_checkbox_binder
"""

# ╔═╡ b51148d5-cce6-4310-b7d4-dcbb6d4ac66b
md"""
Should we plot electrons? $plot_electrons_binder
"""

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

# ╔═╡ dec211fb-33a0-4b16-ad5f-74dc010cfd6f
md"""
# Fitting histograms
"""

# ╔═╡ 4a32f313-8ba7-4354-9843-efd86607efb8
md"""
## Specified bin-width histograms
"""

# ╔═╡ 762fe736-070e-4624-b683-e1fcbeb0f1e0
specific_hist_fits = specific_width_histogram_fits(CR_p_gdf_momentum, 0.01)

# ╔═╡ 80df5127-2df2-4b8d-a2f9-b1d234a96e01
md"""
Proton momentum slice to plot (index): $proton_index_binder (min: $(minimum(idx_CR_p_gdf)), max: $(maximum(idx_CR_p_gdf)))
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
log_p_nat_at_slice = proton_log_p_nat[proton_momentum_index];

# ╔═╡ 89bcb29b-0b1c-4e3a-91cb-282c05df2bc5
md"""
Value of proton momentum at current slice: log(*p*/*m*ₚ*c*) = $log_p_nat_at_slice
"""

# ╔═╡ 4051e244-4c84-4983-8cb9-bc7f53daa9f6
let df = CR_p_gdf_momentum[proton_momentum_index], distribs = normal_distrib_protons
    fig = Figure()
    ax = Axis(
        fig[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",# yscale = log10,
        title = "Histogram of protons dN/dp at log p = $log_p_nat_at_slice (mₚc)",
        axis_properties...)

    if do_plot_pf
        log_dNdp = df.log_dNdp_cr_pf |> skipmissing |> collect
        # log_dNdp ./= std(log_dNdp)
        !isempty(log_dNdp) && stephist!(ax, log_dNdp, label = "plasma frame ($(length(log_dNdp)) samples)"; bins, normalization, color = color_pf_p)

        distrib = distribs.pf[proton_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("MLE fit 𝒩 (%.6f, %.6f)", params(distrib)...), color = :indianred)
        end

        curve_fit_distrib = normal_distrib_protons_from_curves.pf[proton_momentum_index]
        if !ismissing(curve_fit_distrib)
            plot!(ax, curve_fit_distrib, label = @sprintf("Curve fit 𝒩 (%.6f, %.6f)", params(distrib)...), color = :orange)
        end
    end

    if do_plot_sf
        log_dNdp = df.log_dNdp_cr_sf |> skipmissing |> collect
        !isempty(log_dNdp) && hist!(ax, log_dNdp, label = "shock frame"; bins, normalization, color = color_sf_p)
        distrib = distribs.sf[proton_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", params(distrib)...), color = color_sf_p)
        end
    end
    if do_plot_ISM
        log_dNdp = df.log_dNdp_cr_ISM |> skipmissing |> collect
        !isempty(log_dNdp) && hist!(ax, log_dNdp, label = "ISM frame"; bins, normalization, color = color_ISM_p)
        distrib = distribs.ISM[proton_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", params(distrib)...), color = color_ISM_p)
        end
    end

    try
        leg = axislegend(ax, framevisible = false, position = :rt)
        # leg.tellheight = true
    catch e
        # axislegend has no plots to work with, because the current index doesn't have any samples. stop it complaining.
    end
    fig
end

# ╔═╡ 0c230911-62b3-4133-9f17-758bfeb627a2
let
    # bins = 200
    fig = Figure()
    ax = Axis(
        fig[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Truncated Histogram of protons dN/dp at log p = $log_p_nat_at_slice (mₚc)",
        axis_properties...)

    log_dNdp = log_dNdp_cur_trunc
    isempty(log_dNdp) && error("Not the correct momentum slice")
    x, y = get_hist_curve(log_dNdp; nbins=bins)
    # lines!(ax, x, y, label = "bin-centered \"histogram\"", linewidth = 0.5)
    stephist!(ax, log_dNdp, label = "data"; bins, normalization, color = :teal, linewidth = 0.5)
    distrib = fitdistribution(Normal, allowmissing(log_dNdp_cur_trunc))
    ismissing(distrib) || plot!(ax, x, distrib, label = @sprintf("MLE fit 𝒩 (%.2f, %.2f)", params(distrib)...), color = :indianred, linewidth = 1)
    custom_dist = normal_distrib_protons_from_curves.pf[proton_momentum_index]
    plot!(ax, x, custom_dist, label = @sprintf("curve fit 𝒩 (%.2f, %.2f)", params(custom_dist)...), color = :orange, linewidth = 1)

    try
        axislegend(ax, framevisible = false, position = :rt)
    catch e
        # axislegend has no plots to work with, because the current index doesn't have any samples. stop it complaining.
    end
    fig
end

# ╔═╡ 4786bdb6-a387-4333-b9d1-c672dc041910
let
    # bins = 200
    fig = Figure()
    ax = Axis(
        fig[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Histogram (left-tail) of protons dN/dp at log p = $log_p_nat_at_slice (mₚc)",
        axis_properties...)

    log_dNdp = log_dNdp_cur_lowtail
    isempty(log_dNdp) && error("Not the correct momentum slice")
    x, y = get_hist_curve(log_dNdp; nbins=bins)
    # lines!(ax, x, y, label = "bin-centered \"histogram\"", linewidth = 0.5)
    stephist!(ax, log_dNdp, label = "data"; bins, normalization, color = :teal, linewidth = 0.5)
    # lines!(ax, x, fitted_dist_MLE, label = "MLE fit dist", linewidth = 0.5)
    # lines!(ax, x, fitted_dist_curve, label = "curve-fit dist", linewidth = 0.5)
    axislegend(ax, framevisible = false, position = :lt)
    fig
end

# ╔═╡ f7484fdb-37a6-4300-a08d-0e552bc4ef49
let
    # bins = 200
    fig = Figure()
    ax = Axis(
        fig[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Histogram (right-tail) of protons dN/dp at log p = $log_p_nat_at_slice (mₚc)",
        axis_properties...)

    log_dNdp = log_dNdp_cur_hightail
    isempty(log_dNdp) && error("Not the correct momentum slice")
    x, y = get_hist_curve(log_dNdp; nbins=bins)
    # lines!(ax, x, y, label = "bin-centered \"histogram\"", linewidth = 0.5)
    stephist!(ax, log_dNdp, label = "data"; bins, normalization, color = :teal, linewidth = 0.5)
    # distrib = fitdistribution(Normal, allowmissing(log_dNdp_cur_trunc))
    # ismissing(distrib) || plot!(ax, x, distrib, label = @sprintf("MLE fit 𝒩 (%.2f, %.2f)", params(distrib)...), color = :indianred, linewidth = 1)
    # lines!(ax, x, fitted_dist_MLE, label = "MLE fit dist", linewidth = 0.5)
    # lines!(ax, x, fitted_dist_curve, label = "curve-fit dist", linewidth = 0.5)

    axislegend(ax, framevisible = false, position = :rt)
    fig
end

# ╔═╡ 5e67c332-b2cb-45d6-8dff-eeea5acfb779
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Histogram of protons dN/dp at log p = $log_p_nat_at_slice (mₚc)",
        axis_properties...)

    plot!(ax, specific_hist_fits[proton_momentum_index], label = "data"; color = :teal)

    try
        axislegend(ax, framevisible = false, position = :rt)
    catch e
        # axislegend has no plots to work with, because the current index doesn't have any samples. stop it complaining.
    end
    fig
end

# ╔═╡ 589661b1-6a64-4db5-ac40-c1565c29c3cc
electron_log_p_nat = keys(CR_e_gdf_momentum) .|> values .|> first;

# ╔═╡ f186773a-cbd7-47a7-81e4-e0b3e71ebcc5
elec_p_idx = sortperm(electron_log_p_nat)

# ╔═╡ b4df1258-7625-4a7d-b66f-07f4c5e9ba41
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        title = "Sample size vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "# of samples",
    )

    if do_plot_pf
        lengths = gdf_sample_stats(length, CR_p_gdf_momentum; column = :log_dNdp_cr_pf)
        lines!(ax, proton_log_p_nat, lengths, color = color_pf_p, label = "protons, plasma frame")

        if do_plot_electrons
            lengths = gdf_sample_stats(length, CR_e_gdf_momentum; column = :log_dNdp_cr_pf)
            lines!(ax, electron_log_p_nat[elec_p_idx], lengths[elec_p_idx], color = color_pf_e, label = "electrons, plasma frame")
        end
    end
    if do_plot_sf
        lengths = gdf_sample_stats(length, CR_p_gdf_momentum; column = :log_dNdp_cr_sf)
        lines!(ax, proton_log_p_nat, lengths, color = color_sf_p, label = "protons, shock frame")

        if do_plot_electrons
            lengths = gdf_sample_stats(length, CR_e_gdf_momentum; column = :log_dNdp_cr_sf)
            lines!(ax, electron_log_p_nat[elec_p_idx], lengths[elec_p_idx], color = color_sf_e, label = "electrons, shock frame")
        end
    end
    if do_plot_ISM
        lengths = gdf_sample_stats(length, CR_p_gdf_momentum; column = :log_dNdp_cr_ISM)
        lines!(ax, proton_log_p_nat, lengths, color = color_ISM_p, label = "protons, ISM frame")

        if do_plot_electrons
            lengths = gdf_sample_stats(length, CR_e_gdf_momentum; column = :log_dNdp_cr_ISM)
            lines!(ax, electron_log_p_nat[elec_p_idx], lengths[elec_p_idx], color = color_ISM_e, label = "electrons, ISM frame")
        end
    end
    axislegend(ax, framevisible = false, position = :cc)

    fig
end

# ╔═╡ 91bba2da-c925-4123-bb8a-c1f9be8619e9
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        title = "Sample mean vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "⟨log dN/dp⟩",
    )

    if do_plot_pf
        means = gdf_sample_stats(mean, CR_p_gdf_momentum; column = :log_dNdp_cr_pf)
        lines!(ax, proton_log_p_nat, means, color = color_pf_p, label = "protons, plasma frame")

        if do_plot_electrons
            means = gdf_sample_stats(mean, CR_e_gdf_momentum; column = :log_dNdp_cr_pf)
            lines!(ax, electron_log_p_nat, means, color = color_pf_e, label = "electrons, plasma frame")
        end
    end
    if do_plot_sf
        means = gdf_sample_stats(mean, CR_p_gdf_momentum; column = :log_dNdp_cr_sf)
        lines!(ax, proton_log_p_nat, means, color = color_sf_p, label = "protons, shock frame")

        if do_plot_electrons
            means = gdf_sample_stats(mean, CR_e_gdf_momentum; column = :log_dNdp_cr_sf)
            lines!(ax, electron_log_p_nat, means, color = color_sf_e, label = "electrons, shock frame")
        end
    end
    if do_plot_ISM
        means = gdf_sample_stats(mean, CR_p_gdf_momentum; column = :log_dNdp_cr_ISM)
        lines!(ax, proton_log_p_nat, means, color = color_ISM_p, label = "protons, ISM frame")

        if do_plot_electrons
            means = gdf_sample_stats(mean, CR_e_gdf_momentum; column = :log_dNdp_cr_ISM)
            lines!(ax, electron_log_p_nat, means, color = color_ISM_e, label = "electrons, ISM frame")
        end
    end
    axislegend(ax, framevisible = false)

    fig
end

# ╔═╡ 3f36fc06-3799-41f3-971b-d13d43e5fc20
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        title = "Sample mean vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "⟨log dN/dp⟩",
    )

    if do_plot_pf
        means = CR_p_mean_log_dNdp.pf
        std_devs = CR_p_std_log_dNdp.pf
        lines!(ax, proton_log_p_nat, means, color = color_pf_p, label = "protons, plasma frame")
        band!(ax, proton_log_p_nat, means+std_devs, means-std_devs, alpha = 0.4, color = color_pf_p, label = "protons, plasma frame")

        if do_plot_electrons
            means = CR_e_mean_log_dNdp.pf
            std_devs = CR_e_std_log_dNdp.pf
            lines!(ax, electron_log_p_nat, means, color = color_pf_e, label = "electrons, plasma frame")
            band!(ax, electron_log_p_nat, means+std_devs, means-std_devs, alpha = 0.4, color = color_pf_e, label = "electrons, plasma frame")
        end
    end
    if do_plot_sf
        means = CR_p_mean_log_dNdp.sf
        std_devs = CR_p_std_log_dNdp.sf
        lines!(ax, proton_log_p_nat, means, color = color_sf_p, label = "protons, shock frame")
        band!(ax, proton_log_p_nat, means+std_devs, means-std_devs, alpha = 0.4, color = color_sf_p, label = "protons, shock frame")

        if do_plot_electrons
            means = CR_e_mean_log_dNdp.sf
            std_devs = CR_e_std_log_dNdp.sf
            lines!(ax, electron_log_p_nat, means, color = color_sf_e, label = "electrons, shock frame")
            band!(ax, electron_log_p_nat, means+std_devs, means-std_devs, alpha = 0.4, color = color_sf_e, label = "electrons, shock frame")
        end
    end
    if do_plot_ISM
        means = CR_p_mean_log_dNdp.ISM
        std_devs = CR_p_std_log_dNdp.ISM
        lines!(ax, proton_log_p_nat, means, color = color_ISM_p, label = "protons, ISM frame")
        band!(ax, proton_log_p_nat, means+std_devs, means-std_devs, alpha = 0.4, color = color_ISM_p, label = "protons, ISM frame")

        if do_plot_electrons
            means = CR_e_mean_log_dNdp.ISM
            std_devs = CR_e_std_log_dNdp.ISM
            lines!(ax, electron_log_p_nat, means, color = color_ISM_e, label = "electrons, ISM frame")
            band!(ax, electron_log_p_nat, means+std_devs, means-std_devs, alpha = 0.4, color = color_ISM_e, label = "electrons, ISM frame")
        end
    end
    axislegend(ax, framevisible = false, merge = true)

    fig
end

# ╔═╡ 930a7033-e01a-434a-a88a-1bd901dc6bdc
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        title = "Curve-fit σ vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "σ",
        # yscale = log10,
    )
    markersize = 4
    σ_getter = passmissing(d -> d.σ)
    sigmas = σ_getter.(normal_distrib_protons_from_curves.pf)
    scatterlines!(ax, proton_log_p_nat, sigmas, color = color_pf_p, label = "protons, plasma frame"; markersize)
    if do_plot_electrons
        sigmas = σ_getter.(normal_distrib_electrons_from_curves.pf)
        scatterlines!(ax, electron_log_p_nat, sigmas, color = color_pf_e, label = "electrons, plasma frame"; markersize)
    end

    axislegend(ax, position = :lt, framevisible = false)

    fig
end

# ╔═╡ b6ce51e5-b4ff-49eb-83db-ecf3e8a081ac
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        title = "Sample standard deviation vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "σ",
        # yscale = log10,
    )
    markersize = 4

    scatterlines!(ax, proton_log_p_nat, CR_p_std_log_dNdp.pf;
                  color = color_pf_p, label = "protons, plasma frame", markersize)

    if do_plot_electrons
        scatterlines!(ax, electron_log_p_nat, CR_e_std_log_dNdp.pf;
                      color = color_pf_e, label = "electrons, plasma frame", markersize)
    end

    axislegend(ax, position = :lt, framevisible = false)

    fig
end

# ╔═╡ adf24143-4be1-46c7-a63a-fe4dd490791d
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        title = "Sample skewness vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "γ",
        #yscale = log10,
    )

    scatterlines!(ax, proton_log_p_nat, gdf_sample_stats(skewness, CR_p_gdf_momentum; column = :log_dNdp_cr_pf);
                  color = color_pf_p, label = "protons, plasma frame", markersize)
    if do_plot_electrons
        scatterlines!(ax, electron_log_p_nat, gdf_sample_stats(skewness, CR_e_gdf_momentum; column = :log_dNdp_cr_pf);
                      color = color_pf_e, label = "electrons, plasma frame", markersize)
    end

    axislegend(ax, position = :lb, framevisible = false)

    fig
end

# ╔═╡ 67533f87-b016-45fe-b582-53c3c225c056
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        title = "Sample kurtosis vs momentum slice",
        axis_properties...,
        xlabel = "log p (nat)", ylabel = "Kurtosis",
        #yscale = log10,
    )

    scatterlines!(ax, proton_log_p_nat, gdf_sample_stats(kurtosis, CR_p_gdf_momentum; column = :log_dNdp_cr_pf);
                  color = color_pf_p, label = "protons, plasma frame", markersize)
    if do_plot_electrons
        scatterlines!(ax, electron_log_p_nat, gdf_sample_stats(kurtosis, CR_e_gdf_momentum; column = :log_dNdp_cr_pf);
                      color = color_pf_e, label = "electrons, plasma frame", markersize)
    end

    axislegend(ax, position = :lt, framevisible = false)

    fig
end

# ╔═╡ 6c16fc5a-7113-4b6e-abf2-de1275cceda5
log_p_nat_at_slice_e = electron_log_p_nat[electron_momentum_index];

# ╔═╡ c9b9969c-2c7f-436e-b5a1-603138a4e196
md"""
Value of electron momentum at current slice: log(*p*/*m*ₚ*c*) = $log_p_nat_at_slice_e
"""

# ╔═╡ 88822f52-aab8-4931-9091-1909da6c604b
let df = CR_e_gdf_momentum[electron_momentum_index], distribs = normal_distrib_electrons
    fig = Figure()
    ax = Axis(
        fig[1,1];
        xlabel = "log(dN/dp)", ylabel = "pdf",
        title = "Histogram of electrons dN/dp at log p = $log_p_nat_at_slice_e (mₚc)",
        axis_properties...)

    if do_plot_pf
        log_dNdp = df.log_dNdp_cr_pf |> skipmissing |> collect
        !isempty(log_dNdp) && stephist!(ax, log_dNdp, label = "plasma frame ($(length(log_dNdp)) samples)"; bins, normalization, color = color_pf_e)
        distrib = distribs.pf[electron_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("MLE fit 𝒩 (%.2f, %.2f)", params(distrib)...), color = :lightgoldenrod4)
        end
        distrib = normal_distrib_electrons_from_curves.pf[electron_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("Curve fit 𝒩 (%.2f, %.2f)", params(distrib)...), color = :green)
        end
    end

    if do_plot_sf
        log_dNdp = df.log_dNdp_cr_sf |> skipmissing |> collect
        !isempty(log_dNdp) && stephist!(ax, log_dNdp, label = "shock frame"; bins, normalization, color = color_sf_e)
        distrib = distribs.sf[electron_momentum_index]
        if !ismissing(distrib)
            μ, σ = params(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", μ, σ), color = color_sf_e)
        end
    end
    if do_plot_ISM
        log_dNdp = df.log_dNdp_cr_ISM |> skipmissing |> collect
        !isempty(log_dNdp) && stephist!(ax, log_dNdp, label = "ISM frame"; bins, normalization, color = color_ISM_e)
        distrib = distribs.ISM[electron_momentum_index]
        if !ismissing(distrib)
            plot!(ax, distrib, label = @sprintf("𝒩 (%.2f, %.2f)", params(distrib)...), color = color_ISM_e)
        end
    end

    try
        axislegend(ax, framevisible = false, position = :lt)
    catch e
        # axislegend has no plots to work with, because the current index doesn't have any samples. stop it complaining.
    end
    fig
end

# ╔═╡ 7534104f-885d-48c5-8ae0-ddae56fcd86d
let
    fig = Figure()
    ax = Axis(fig[1,1];
              yscale = log10,
              axis_properties...,
              title = "Distribution agreement curve",
              xlabel = "log p (nat)", ylabel = "Bhattacharya distance")
    scatterlines!(ax, proton_log_p_nat, proton_distances; label = "protons, plasma frame", markersize)
    scatterlines!(ax, electron_log_p_nat, electron_distances; label = "electrons, plasma frame", markersize)
    axislegend(ax, position = :ct)
    fig
end

# ╔═╡ e7a26d10-0e00-444d-a8f9-27874a8f821e
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
        title = "Root-Sum-of-Squared-Errors vs momentum slice",
        axis_properties...,
        xminorticksvisible = true, yminorticksvisible = true,
        xlabel = "log p (nat)",
        yscale = p_val_yscale,
    )

    # scatterlines!(ax, proton_log_p_nat, sse_scores_p, color = color_pf_p, label = "protons, plasma frame (MLE)"; markersize)
    scatterlines!(ax, proton_log_p_nat, sse_scores_p_curve, label = "protons, plasma frame (curve)"; markersize, color = :darkgreen)
    if do_plot_electrons
        # scatterlines!(ax, electron_log_p_nat, sse_scores_e, color = color_pf_e, label = "electrons, plasma frame (MLE)"; markersize)
        scatterlines!(ax, electron_log_p_nat, sse_scores_e_curve, color = :orange, label = "electrons, plasma frame (curve)"; markersize)
    end

    axislegend(ax, position = plot_p_values_in_logscale ? :ct : :cb, framevisible = false)

    fig
end

# ╔═╡ cee91c99-adc0-4185-a7c3-e2164b95a003
let
    fig = Figure()
    ax = Axis(
        fig[1,1];
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
        fig[1,1];
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
        fig[1,1];
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
# ╟─a36ea9cf-176f-40bd-8577-cc2ea8db64af
# ╟─985a2460-3fbc-4935-af59-2e734786c973
# ╟─d85427f4-86ed-4c04-980a-a4152b5875e8
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
# ╠═4553a97d-6b78-4268-90de-d8bee348d3d4
# ╟─3bf64608-0fa2-4fcb-9782-fd7a8de47bda
# ╠═70717f68-e97b-401e-bcf7-0684ade30b07
# ╠═c09cf00d-0a07-4159-9009-45afdb8343fb
# ╠═a8fdec46-2d8d-4a06-800a-7e03c64468dd
# ╠═a9b9b9aa-b850-49de-8412-2930e7004a36
# ╠═734626db-7646-4a0c-9303-2aec6e371159
# ╟─e85bc3ef-97a4-4d0e-937e-3ca290a28086
# ╟─88527827-9710-4c63-964f-691aa8909d1f
# ╟─cd684b62-1aad-404a-87b4-e7b406c0c989
# ╟─54ec7213-3e48-423b-9208-befc6583b908
# ╠═f186773a-cbd7-47a7-81e4-e0b3e71ebcc5
# ╟─b4df1258-7625-4a7d-b66f-07f4c5e9ba41
# ╟─932c2a77-0198-4df4-a4bd-30d0bda93946
# ╟─91bba2da-c925-4123-bb8a-c1f9be8619e9
# ╟─7da22399-718c-4f16-ba02-7ae27773942b
# ╟─3f36fc06-3799-41f3-971b-d13d43e5fc20
# ╟─5767b9ac-64c2-4d2f-ad42-961184c7edc7
# ╟─930a7033-e01a-434a-a88a-1bd901dc6bdc
# ╟─b6ce51e5-b4ff-49eb-83db-ecf3e8a081ac
# ╟─7495e7e9-3d50-4401-baef-d2e3c11e6b46
# ╟─adf24143-4be1-46c7-a63a-fe4dd490791d
# ╟─bb3a74ce-78ee-487e-a413-4c0e035e8818
# ╟─67533f87-b016-45fe-b582-53c3c225c056
# ╟─ce8b1307-dc78-463b-9f41-04fe5dded525
# ╠═71404de8-f8b2-4d26-b7d7-41064cae1447
# ╠═6c16fc5a-7113-4b6e-abf2-de1275cceda5
# ╟─89bcb29b-0b1c-4e3a-91cb-282c05df2bc5
# ╟─c9b9969c-2c7f-436e-b5a1-603138a4e196
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
# ╠═febdc8a1-00bb-47a7-83d2-6cccef5190f5
# ╟─35710ad9-f2e4-487b-be19-c29500633726
# ╟─4051e244-4c84-4983-8cb9-bc7f53daa9f6
# ╟─7be1e6da-0eb9-45e5-a4f9-bb6deedc3def
# ╟─88822f52-aab8-4931-9091-1909da6c604b
# ╟─9ea7a3a4-987d-416d-88d1-672e3cce23c5
# ╟─b7a96870-784e-4ce0-830d-d245fc16e5f4
# ╟─4ac1798d-ec27-4571-9b2a-44cb432ef0d6
# ╟─a6a63cb1-1a13-4cc2-9730-b78dd3d3aee4
# ╠═54452e38-227e-4d06-ae74-7347aae2c021
# ╠═b0d555b3-5087-4405-8343-ce304d482ca9
# ╠═5dc367ca-2882-4b98-8f29-2b5390426a9b
# ╟─464e92d7-e414-4ddd-a81e-978f271961b2
# ╟─daded5be-4319-4a0a-8491-84a7273a844b
# ╠═5bbd6e99-87e1-401c-a09e-065e2d426370
# ╠═5ab05dc9-3a98-4297-a47b-c4e0111b8c51
# ╟─aebbf1a7-d047-4fa6-aa36-4b9ae8b68127
# ╠═55d8c831-27e6-4914-a836-7a05281e8fb3
# ╠═89f8d7a8-ea2e-4906-9460-da16154b0404
# ╠═29ec59ad-0e22-462a-ab6d-2065a56fc001
# ╠═32f07cd2-f62f-41e0-9211-8ac333bdd98d
# ╟─6cb898b3-98c5-4f3a-8d77-3deef7cf5358
# ╟─2374b968-1172-48db-8ddd-7b4deae7817c
# ╠═c79d7166-e0a5-465e-b78d-6eeee331a99c
# ╠═a2c5161c-7cbd-4314-8fda-59a8e1750da1
# ╠═e5d14c88-fea2-4117-97ab-0aebf3711a5b
# ╠═59444b54-893e-4f4e-b746-97de78417043
# ╠═316b7ed7-4fe8-40dc-8b1a-551a0100c57c
# ╠═07300c09-2361-40e0-a502-b018496184c8
# ╠═afabc297-408f-4643-8296-40be885adafc
# ╠═60bd9873-0246-432d-9e68-bfe2aa0956b2
# ╟─46607c86-7310-4d81-9816-2283d28d1420
# ╟─0c230911-62b3-4133-9f17-758bfeb627a2
# ╟─4786bdb6-a387-4333-b9d1-c672dc041910
# ╟─f7484fdb-37a6-4300-a08d-0e552bc4ef49
# ╟─f3132403-113d-4b30-9fd0-379d28ade3c7
# ╟─8bc20375-2562-4611-a67b-5884aa99b5f0
# ╠═e6b9701d-3d27-4c0c-b0b9-9879527f369c
# ╠═e75ea9c0-59ca-4097-b4f6-6a3af04dc308
# ╟─639ab710-c411-4831-9e1d-d7fba723b7bc
# ╟─3860d0cf-20f4-4256-9286-8757afc38ef9
# ╠═f2fe3844-2be8-4da6-9656-40312304556b
# ╠═97291776-74f0-428a-ab4f-3c498b630000
# ╟─452c9b2f-7138-4310-b0c6-df2be7ab8c76
# ╠═b238afe1-3d1f-4e15-bc49-1b015a39c02c
# ╠═78a22648-c76a-4b5c-b552-9be000a60109
# ╠═7534104f-885d-48c5-8ae0-ddae56fcd86d
# ╟─da107273-c428-4c68-80a9-8f82cb211497
# ╠═49902e99-870d-4d19-afb0-1de612c185df
# ╠═4ac32bef-af81-4f7e-8e97-7eac4dd2bf69
# ╟─f330af91-60a6-46ac-bdc5-ec49c216fccb
# ╠═fd47dab7-426c-44fc-8038-00e378324e41
# ╠═d77a95bf-2d54-46d5-81ca-b671ee1db695
# ╟─04dad413-0dc0-4ceb-81c2-e208ef082f38
# ╟─94a91acd-a878-4c3c-9716-8bed60bf8c6c
# ╠═222df0cb-0760-48a2-902e-91d32e451a11
# ╠═70f93b37-a977-4dce-9fdd-a0497603a864
# ╠═c5d56d28-15d5-464f-9055-7bbfc9826e72
# ╠═cbea4ff4-b132-4abb-97c6-e406a339ced6
# ╠═85feafa4-a572-40ca-9975-fb0d3d5309f7
# ╟─79dc57bb-d66d-4608-a775-9dfc58af1995
# ╟─b51148d5-cce6-4310-b7d4-dcbb6d4ac66b
# ╟─e7a26d10-0e00-444d-a8f9-27874a8f821e
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
# ╟─dec211fb-33a0-4b16-ad5f-74dc010cfd6f
# ╟─4a32f313-8ba7-4354-9843-efd86607efb8
# ╠═762fe736-070e-4624-b683-e1fcbeb0f1e0
# ╟─80df5127-2df2-4b8d-a2f9-b1d234a96e01
# ╠═5e67c332-b2cb-45d6-8dff-eeea5acfb779
# ╟─8d03de5e-d344-4efd-b9af-dd5391028780
# ╠═377aaf8f-b909-4c42-bc77-912fd300c300
# ╟─32edc221-e586-4510-9427-977b22f62f6c
# ╠═e8406a6a-ecc2-49d2-b67a-503b4ef5764b
# ╠═589661b1-6a64-4db5-ac40-c1565c29c3cc
