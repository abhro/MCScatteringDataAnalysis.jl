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
    should_title = args["title"]

    CR_p_gdf_momentum_filename = joinpath(dir, "dNdp-CR-protons-momentum-split.jld2")

    @info("Starting plots of protons")
    species_plots(CR_p_gdf_momentum_filename; species_name = "protons", bins, outdir, should_title)
    @info("Finished plots of protons")

    GC.gc()

    CR_e_gdf_momentum_filename = joinpath(dir, "dNdp-CR-electrons-momentum-split.jld2")

    @info("Starting plots of electrons")
    species_plots(CR_e_gdf_momentum_filename; species_name = "electrons", bins, outdir, should_title)
    @info("Finished plots of electrons")

    return
end

function species_plots(CR_gdf_filename; species_name, bins, outdir, should_title)
    (log_p_nat, mle_fit_dists, hist_fit_dists) = dist_fits(CR_gdf_filename, bins)
    @info("Fitted $species_name data to distributions")
    @debug("Received fits", mle_fit_dists, hist_fit_dists)
    fig = Figure()
    ax = Axis(fig[1, 1]; axis_properties...)
    if should_title
        ax.title = plot_title
    end
    hlines!(ax, 2.0e-3, linewidth = 0.8, linestyle = :dash)
    for (bin_count, dist) in zip(bins, hist_fit_dists)
        @debug("Processing bin count=$bin_count", dist)
        distances = bcdistances(mle_fit_dists.pf, dist.pf)
        scatterlines!(ax, log_p_nat, distances; label = "bins = $bin_count", markersize)
    end
    axislegend(ax, position = :ct, framevisible = false)
    plot_filename = "bc-distances-$species_name-$(join(bins, ","))-bins.svg"
    save(joinpath(outdir, plot_filename), fig)
    return
end

"""
Read a file, for each momentum slice, make an MLE fit and a bunch of histogram
fits based on number of bins provided
"""
function dist_fits(filename, bins)
    df = load_object(filename)
    log_p = keys(df) .|> values .|> first
    # MLE fitted distribution
    mle = fitdistributions(v -> fitdistribution(Normal, v), df)
    # distributions fitted through least squares on histogram curve
    hist_dists = []
    for bin_count in bins
        push!(
            hist_dists,
            fitdistributions(
                v -> fit_dist_to_histogram(Normal, v; nbins = bin_count),
                df,
            )
        )
    end
    return (; log_p, mle, hist_dists)
end


const plot_title = "Distribution agreement curve"
const axis_properties = (;
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorticksvisible = true,
    yminorticksvisible = true,
    xlabel = "log p (nat)",
    ylabel = "Bhattacharya distance",
    yscale = log10,
)

const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors()

"""
Create and return an ArgParse.jl arg table
"""
function get_parser()
    s = ArgParseSettings(
        help_width = 120
    )

    add_arg_table!(
        s,
        "data-dir",
        Dict(:help => "Where to look for the .jld2 files", :required => true),
        "out-dir",
        Dict(:help => "Where to save the plots", :required => true),

        "--bins",
        Dict(
            :help => "Number of bins to use for the histogram",
            :default => [90],
            :arg_type => Int64,
            :nargs => '+',
            :metavar => "N",
        ),
        "--column",
        Dict(:help => "Which column of the dataframe to bin", :default => :log_dNdp_cr_pf),

        "--log-scale",
        Dict(:help => "Whether the y-axis should be in log-scale", :action => :store_true),

        "--title",
        Dict(:help => "Whether we should add a title to the figure", :action => :store_true),
    )
    return s
end
