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

# ╔═╡ e5e0e4e2-2df1-4536-9cc5-bdcec6fc13de
using DrWatson: datadir

# ╔═╡ 7899ae97-fbc2-43e5-ac77-c6d725f0371e
using JLD2, DataFrames

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

# ╔═╡ f0e77bbd-e420-49f1-9b40-f9d994888b93
md"""
# Plot flux statistics for each momentum slice
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

# ╔═╡ efed6856-264e-4917-8d01-aadfa4aaf9ea
proton_log_p_nat = keys(CR_p_gdf_momentum) .|> values .|> first;

# ╔═╡ 72e45c0b-d1e4-4244-8ff8-0d1013a57d71
electron_log_p_nat = keys(CR_e_gdf_momentum) .|> values .|> first;

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
    xlabel = "log p (nat)",
)

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

# ╔═╡ 4553a97d-6b78-4268-90de-d8bee348d3d4
plot_electrons_binder = @bind do_plot_electrons CheckBox(default = true);

# ╔═╡ 3bf64608-0fa2-4fcb-9782-fd7a8de47bda
md"""
## Sample statistics
"""

# ╔═╡ ad0e0d21-707c-4c3f-bc5f-3e3718421f13
function span(v::AbstractArray{T}) where {T}
    if isempty(v)
        return zero(T)
    end
    vmin, vmax = extrema(v)
    return vmax - vmin
end

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

# ╔═╡ 8b453952-b190-45c3-9cfa-5a5df29bbe0f
CR_p_std_scaled_log_dNdp = (;
    pf = gdf_sample_stats(v -> std(v) / span(v), CR_p_gdf_momentum; column = :log_dNdp_cr_pf),
    sf = gdf_sample_stats(v -> std(v) / span(v), CR_p_gdf_momentum; column = :log_dNdp_cr_sf),
    ISM = gdf_sample_stats(v -> std(v) / span(v), CR_p_gdf_momentum; column = :log_dNdp_cr_ISM),
);

# ╔═╡ 7a17522e-6ab4-48bc-9902-81b760ec01b1
CR_e_std_scaled_log_dNdp = (;
    pf = gdf_sample_stats(v -> std(v) / span(v), CR_e_gdf_momentum; column = :log_dNdp_cr_pf),
    sf = gdf_sample_stats(v -> std(v) / span(v), CR_e_gdf_momentum; column = :log_dNdp_cr_sf),
    ISM = gdf_sample_stats(v -> std(v) / span(v), CR_e_gdf_momentum; column = :log_dNdp_cr_ISM),
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
## Sample size
"""

# ╔═╡ f186773a-cbd7-47a7-81e4-e0b3e71ebcc5
elec_p_idx = sortperm(electron_log_p_nat)

# ╔═╡ b4df1258-7625-4a7d-b66f-07f4c5e9ba41
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Sample size vs momentum slice",
        axis_properties...,
        ylabel = "# of samples",
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

# ╔═╡ 932c2a77-0198-4df4-a4bd-30d0bda93946
md"""
## Means
"""

# ╔═╡ 91bba2da-c925-4123-bb8a-c1f9be8619e9
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Sample mean vs momentum slice",
        axis_properties...,
        ylabel = "⟨log dN/dp⟩",
    )

    if do_plot_pf
        means = CR_p_mean_log_dNdp.pf
        lines!(ax, proton_log_p_nat, means, color = color_pf_p, label = "protons, plasma frame")

        if do_plot_electrons
            means = CR_e_mean_log_dNdp.pf
            lines!(ax, electron_log_p_nat, means, color = color_pf_e, label = "electrons, plasma frame")
        end
    end
    if do_plot_sf
        means = CR_p_mean_log_dNdp.sf
        lines!(ax, proton_log_p_nat, means, color = color_sf_p, label = "protons, shock frame")

        if do_plot_electrons
            means = CR_e_mean_log_dNdp.sf
            lines!(ax, electron_log_p_nat, means, color = color_sf_e, label = "electrons, shock frame")
        end
    end
    if do_plot_ISM
        means = CR_p_mean_log_dNdp.ISM
        lines!(ax, proton_log_p_nat, means, color = color_ISM_p, label = "protons, ISM frame")

        if do_plot_electrons
            means = CR_e_mean_log_dNdp.ISM
            lines!(ax, electron_log_p_nat, means, color = color_ISM_e, label = "electrons, ISM frame")
        end
    end
    axislegend(ax, framevisible = false)

    fig
end

# ╔═╡ 7da22399-718c-4f16-ba02-7ae27773942b
md"""
### Mean with uncertainty envelope
"""

# ╔═╡ 3f36fc06-3799-41f3-971b-d13d43e5fc20
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Sample mean vs momentum slice",
        axis_properties...,
        ylabel = "⟨log dN/dp⟩",
    )

    if do_plot_pf
        means = CR_p_mean_log_dNdp.pf
        std_devs = CR_p_std_log_dNdp.pf
        lines!(ax, proton_log_p_nat, means, color = color_pf_p, label = "protons, plasma frame")
        band!(ax, proton_log_p_nat, means + std_devs, means - std_devs, alpha = 0.4, color = color_pf_p, label = "protons, plasma frame")

        if do_plot_electrons
            means = CR_e_mean_log_dNdp.pf
            std_devs = CR_e_std_log_dNdp.pf
            lines!(ax, electron_log_p_nat, means, color = color_pf_e, label = "electrons, plasma frame")
            band!(ax, electron_log_p_nat, means + std_devs, means - std_devs, alpha = 0.4, color = color_pf_e, label = "electrons, plasma frame")
        end
    end
    if do_plot_sf
        means = CR_p_mean_log_dNdp.sf
        std_devs = CR_p_std_log_dNdp.sf
        lines!(ax, proton_log_p_nat, means, color = color_sf_p, label = "protons, shock frame")
        band!(ax, proton_log_p_nat, means + std_devs, means - std_devs, alpha = 0.4, color = color_sf_p, label = "protons, shock frame")

        if do_plot_electrons
            means = CR_e_mean_log_dNdp.sf
            std_devs = CR_e_std_log_dNdp.sf
            lines!(ax, electron_log_p_nat, means, color = color_sf_e, label = "electrons, shock frame")
            band!(ax, electron_log_p_nat, means + std_devs, means - std_devs, alpha = 0.4, color = color_sf_e, label = "electrons, shock frame")
        end
    end
    if do_plot_ISM
        means = CR_p_mean_log_dNdp.ISM
        std_devs = CR_p_std_log_dNdp.ISM
        lines!(ax, proton_log_p_nat, means, color = color_ISM_p, label = "protons, ISM frame")
        band!(ax, proton_log_p_nat, means + std_devs, means - std_devs, alpha = 0.4, color = color_ISM_p, label = "protons, ISM frame")

        if do_plot_electrons
            means = CR_e_mean_log_dNdp.ISM
            std_devs = CR_e_std_log_dNdp.ISM
            lines!(ax, electron_log_p_nat, means, color = color_ISM_e, label = "electrons, ISM frame")
            band!(ax, electron_log_p_nat, means + std_devs, means - std_devs, alpha = 0.4, color = color_ISM_e, label = "electrons, ISM frame")
        end
    end
    axislegend(ax, framevisible = false, merge = true)

    fig
end

# ╔═╡ 5767b9ac-64c2-4d2f-ad42-961184c7edc7
md"""
## Standard deviations
"""

# ╔═╡ b6ce51e5-b4ff-49eb-83db-ecf3e8a081ac
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Sample standard deviation vs momentum slice",
        axis_properties...,
        ylabel = "σ",
        # yscale = log10,
    )
    markersize = 4

    if do_plot_pf
        std_devs = CR_p_std_log_dNdp.pf
        scatterlines!(ax, proton_log_p_nat, std_devs; color = color_pf_p, label = "protons, plasma frame", markersize)

        if do_plot_electrons
            std_devs = CR_e_std_log_dNdp.pf
            scatterlines!(ax, electron_log_p_nat, std_devs; color = color_pf_e, label = "electrons, plasma frame", markersize)
        end
    end
    if do_plot_sf
        std_devs = CR_p_std_log_dNdp.sf
        scatterlines!(ax, proton_log_p_nat, std_devs; color = color_sf_p, label = "protons, shock frame", markersize)

        if do_plot_electrons
            std_devs = CR_e_std_log_dNdp.sf
            scatterlines!(ax, electron_log_p_nat, std_devs; color = color_sf_e, label = "electrons, shock frame", markersize)
        end
    end
    if do_plot_ISM
        std_devs = CR_p_std_log_dNdp.ISM
        scatterlines!(ax, proton_log_p_nat, std_devs; color = color_ISM_p, label = "protons, ISM frame", markersize)

        if do_plot_electrons
            std_devs = CR_e_std_log_dNdp.ISM
            scatterlines!(ax, electron_log_p_nat, std_devs; color = color_ISM_e, label = "electrons, ISM frame", markersize)
        end
    end
    axislegend(ax, position = :lt, framevisible = false)

    fig
end

# ╔═╡ 6266b083-269c-442b-8682-2587c944a276
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Sample standard deviation over data range vs momentum slice",
        axis_properties...,
        ylabel = "σ/span(v)",
        yscale = log10,
    )
    markersize = 4

    if do_plot_pf
        scatterlines!(ax, proton_log_p_nat, CR_p_std_scaled_log_dNdp.pf; color = color_pf_p, label = "protons, plasma frame", markersize)

        if do_plot_electrons
            scatterlines!(ax, electron_log_p_nat, CR_e_std_scaled_log_dNdp.pf; color = color_pf_e, label = "electrons, plasma frame", markersize)
        end
    end
    if do_plot_sf
        scatterlines!(ax, proton_log_p_nat, CR_p_std_scaled_log_dNdp.sf; color = color_sf_p, label = "protons, shock frame", markersize)

        if do_plot_electrons
            scatterlines!(ax, electron_log_p_nat, CR_e_std_scaled_log_dNdp.sf; color = color_sf_e, label = "electrons, shock frame", markersize)
        end
    end
    if do_plot_ISM
        scatterlines!(ax, proton_log_p_nat, CR_p_std_scaled_log_dNdp.ISM; color = color_ISM_p, label = "protons, ISM frame", markersize)

        if do_plot_electrons
            scatterlines!(ax, electron_log_p_nat, CR_e_std_scaled_log_dNdp.ISM; color = color_ISM_e, label = "electrons, ISM frame", markersize)
        end
    end
    axislegend(ax, position = :lt, framevisible = false)

    fig
end

# ╔═╡ 7495e7e9-3d50-4401-baef-d2e3c11e6b46
md"""
## Skewness
"""

# ╔═╡ de5f1493-a018-40cd-8b30-0681f1c61768
md"""
Should we plot electrons? $plot_electrons_binder
"""

# ╔═╡ adf24143-4be1-46c7-a63a-fe4dd490791d
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Sample skewness vs momentum slice",
        axis_properties...,
        ylabel = "γ",
        #yscale = log10,
    )

    scatterlines!(
        ax, proton_log_p_nat, gdf_sample_stats(skewness, CR_p_gdf_momentum; column = :log_dNdp_cr_pf);
        color = color_pf_p, label = "protons, plasma frame", markersize
    )
    if do_plot_electrons
        scatterlines!(
            ax, electron_log_p_nat, gdf_sample_stats(skewness, CR_e_gdf_momentum; column = :log_dNdp_cr_pf);
            color = color_pf_e, label = "electrons, plasma frame", markersize
        )
    end

    axislegend(ax, position = :lb, framevisible = false)

    fig
end

# ╔═╡ bb3a74ce-78ee-487e-a413-4c0e035e8818
md"""
## Kurtosis
"""

# ╔═╡ 67533f87-b016-45fe-b582-53c3c225c056
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Sample kurtosis vs momentum slice",
        axis_properties...,
        ylabel = "Kurtosis",
        #yscale = log10,
    )

    scatterlines!(
        ax, proton_log_p_nat, gdf_sample_stats(kurtosis, CR_p_gdf_momentum; column = :log_dNdp_cr_pf);
        color = color_pf_p, label = "protons, plasma frame", markersize
    )
    if do_plot_electrons
        scatterlines!(
            ax, electron_log_p_nat, gdf_sample_stats(kurtosis, CR_e_gdf_momentum; column = :log_dNdp_cr_pf);
            color = color_pf_e, label = "electrons, plasma frame", markersize
        )
    end

    axislegend(ax, position = :lt, framevisible = false)

    fig
end

# ╔═╡ ecf80697-b786-4b02-9563-f3d082383b76
md"""
Choose which frames to plot:
- Plasma frame: $plot_pf_binder
- Shock frame: $plot_sf_binder
- ISM frame: $plot_ISM_binder
"""

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
# ╠═efed6856-264e-4917-8d01-aadfa4aaf9ea
# ╠═72e45c0b-d1e4-4244-8ff8-0d1013a57d71
# ╟─628130bf-da25-4799-8e5e-3d2db15b1e49
# ╟─ecf233ad-d75e-4aa5-bf7e-ff3e7b1d8755
# ╟─59a22149-3397-4e97-9f7b-5d502aacf293
# ╟─f91132bd-28af-4a6c-9a77-5c5b0ed4a08a
# ╠═50b1a87f-49ff-4d93-aa6e-f042a87b875e
# ╟─105361e9-cafd-4755-bcbd-fdcbcb07b291
# ╟─1b9f507c-1585-4ad1-8090-bdde6de972d6
# ╠═d6513272-7232-43e6-ac88-a58462181041
# ╠═59be6983-6e37-4a70-8929-69176a5f807e
# ╠═60deb76f-3efe-4e0d-b176-9f0169259dca
# ╠═4553a97d-6b78-4268-90de-d8bee348d3d4
# ╟─3bf64608-0fa2-4fcb-9782-fd7a8de47bda
# ╠═ad0e0d21-707c-4c3f-bc5f-3e3718421f13
# ╟─e85bc3ef-97a4-4d0e-937e-3ca290a28086
# ╟─88527827-9710-4c63-964f-691aa8909d1f
# ╟─cd684b62-1aad-404a-87b4-e7b406c0c989
# ╟─54ec7213-3e48-423b-9208-befc6583b908
# ╠═f186773a-cbd7-47a7-81e4-e0b3e71ebcc5
# ╟─b4df1258-7625-4a7d-b66f-07f4c5e9ba41
# ╟─932c2a77-0198-4df4-a4bd-30d0bda93946
# ╠═c09cf00d-0a07-4159-9009-45afdb8343fb
# ╠═a8fdec46-2d8d-4a06-800a-7e03c64468dd
# ╟─91bba2da-c925-4123-bb8a-c1f9be8619e9
# ╟─7da22399-718c-4f16-ba02-7ae27773942b
# ╟─3f36fc06-3799-41f3-971b-d13d43e5fc20
# ╟─5767b9ac-64c2-4d2f-ad42-961184c7edc7
# ╠═a9b9b9aa-b850-49de-8412-2930e7004a36
# ╠═8b453952-b190-45c3-9cfa-5a5df29bbe0f
# ╠═b6ce51e5-b4ff-49eb-83db-ecf3e8a081ac
# ╠═b6ce51e5-b4ff-49eb-83db-ecf3e8a081ac
# ╠═7a17522e-6ab4-48bc-9902-81b760ec01b1
# ╠═734626db-7646-4a0c-9303-2aec6e371159
# ╠═6266b083-269c-442b-8682-2587c944a276
# ╟─7495e7e9-3d50-4401-baef-d2e3c11e6b46
# ╟─de5f1493-a018-40cd-8b30-0681f1c61768
# ╟─adf24143-4be1-46c7-a63a-fe4dd490791d
# ╟─bb3a74ce-78ee-487e-a413-4c0e035e8818
# ╟─67533f87-b016-45fe-b582-53c3c225c056
# ╟─ecf80697-b786-4b02-9563-f3d082383b76
