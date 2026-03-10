using JLD2: load_object
using DataFrames
using LinearAlgebra: normalize
using CairoMakie
using StatsBase
using Distributions: Normal
using ArgParse: ArgParseSettings, add_arg_table!, parse_args
using MCScatteringDataAnalysis: fitdistribution, fitdistributions, fit_dist_to_histogram, bcdistances

const normalization = :pdf
const markersize = 5

"""
Entry point of the program/script.
"""
function (@main)(args)

    # overwrite args local with its parsed version
    args = parse_args(args, get_parser())

    dir = args["data-dir"]
    outdir = args["out-dir"]    # directory where plots will be saved
    @info("Ensuring directory $outdir exists")
    mkpath(outdir)

    column = Symbol(args["column"])
    bins = args["bins"]

    yscale = args["log-scale"] ? log10 : identity

    CR_p_gdf_momentum_filename = joinpath(dir, "dNdp-CR-protons-momentum-split.jld2")

    @info("Starting plots of protons")
    proton_plots(CR_p_gdf_momentum_filename; bins, outdir)
    @info("Finished plots of protons")

    GC.gc()

    CR_e_gdf_momentum_filename = joinpath(dir, "dNdp-CR-electrons-momentum-split.jld2")

    @info("Starting plots of electrons")
    electron_plots(CR_e_gdf_momentum_filename; bins, outdir)
    @info("Finished plots of electrons")

    return
end

function proton_plots(CR_p_gdf_momentum_filename; bins, outdir)
    CR_p_gdf_momentum = load_object(CR_p_gdf_momentum_filename)
    proton_log_p_nat = keys(CR_p_gdf_momentum) .|> values .|> first
    # create histogram for each of the momentum slices
    normal_distrib_protons_from_curves = fitdistributions(v -> fit_dist_to_histogram(Normal, v; nbins=bins), CR_p_gdf_momentum)
    # create MLE fitted distribution for each of the slices
    normal_distrib_protons = fitdistributions(v -> fitdistribution(Normal, v), CR_p_gdf_momentum)
    proton_distances = bcdistances(normal_distrib_protons.pf, normal_distrib_protons_from_curves.pf)
    fig = Figure()
    ax = Axis(fig[1, 1]; axis_properties...)
    scatterlines!(ax, proton_log_p_nat, proton_distances; label = "protons, plasma frame", markersize, color = color_pf_p)
    hlines!(ax, 1.0e-3, linewidth = 0.4, color = color_pf_p, linestyle = :dash)
    axislegend(ax, position = :ct, framevisible = false)
    save(joinpath(outdir, "bc-distances-protons-$bins-bins.svg"), fig)
end

function electron_plots(CR_e_gdf_momentum_filename; bins, outdir)
    CR_e_gdf_momentum = load_object(CR_e_gdf_momentum_filename)
    electron_log_p_nat = keys(CR_e_gdf_momentum) .|> values .|> first
    normal_distrib_electrons = fitdistributions(v -> fitdistribution(Normal, v), CR_e_gdf_momentum)
    normal_distrib_electrons_from_curves = fitdistributions(v -> fit_dist_to_histogram(Normal, v; nbins = bins), CR_e_gdf_momentum)
    electron_distances = bcdistances(normal_distrib_electrons.pf, normal_distrib_electrons_from_curves.pf)
    fig = Figure()
    ax = Axis(fig[1, 1]; axis_properties...)
    scatterlines!(ax, electron_log_p_nat, electron_distances; label = "electrons, plasma frame", markersize, color = color_pf_e)
    hlines!(ax, 1.0e-3, linewidth = 0.4, color = color_pf_e, linestyle = :dash)
    axislegend(ax, position = :ct, framevisible = false)
    save(joinpath(outdir, "bc-distances-electrons-$bins-bins.svg"), fig)
end

const axis_properties = (;
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorticksvisible = true,
    yminorticksvisible = true,
    title = "Distribution agreement curve",
    xlabel = "log p (nat)",
    ylabel = "Bhattacharya distance",
    yscale = log10,
)

const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors()

"""
Create and return an ArgParse.jl arg table
"""
function get_parser()
    s = ArgParseSettings()

    add_arg_table!(
        s,
        "data-dir",
        Dict(:help => "Where to look for the .jld2 files", :required => true),
        "out-dir",
        Dict(:help => "Where to save the plots", :required => true),

        "--bins",
        Dict(
            :help => "Number of bins to use for the histogram",
            :default => 90,
            :arg_type => Int64,
            :metavar => "N",
        ),
        "--column",
        Dict(:help => "Which column of the dataframe to bin", :default => :log_dNdp_cr_pf),

        "--log-scale",
        Dict(:help => "Whether the y-axis should be in log-scale", :action => :store_true),
    )
    return s
end

