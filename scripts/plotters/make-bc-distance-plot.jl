using JLD2: load_object
using DataFrames
using LinearAlgebra: normalize
using CairoMakie
using StatsBase
using Distributions: Normal
using ArgParse: ArgParseSettings, add_arg_table!, parse_args
using MCScatteringDataAnalysis: fitdistribution, fitdistributions, fit_dist_to_histogram, bcdistances, unit_ticks

const normalization = :pdf
const markersize = 5

"""
Entry point of the program/script.
"""
function (@main)(args = [])

    # overwrite args local with its parsed version
    args = parse_args(args, get_parser())

    dir = args["data-dir"]
    outdir = args["out-dir"]    # directory where plots will be saved
    @info("Ensuring directory $outdir exists")
    mkpath(outdir)

    bins = args["bins"]

    yscale = args["log-scale"] ? log10 : identity
    should_title = args["title"]

    Makie.set_theme!(theme_latexfonts())
    Makie.update_theme!(fonts = Attributes(
        regular = "Libertinus Serif",
        bold = "Libertinus Serif Bold",
    ))

    CR_p_gdf_momentum_filename = joinpath(dir, "dNdp-CR-protons-momentum-split.jld2")
    process_species(CR_p_gdf_momentum_filename, "protons"; bins, outdir, should_title)

    GC.gc()

    CR_e_gdf_momentum_filename = joinpath(dir, "dNdp-CR-electrons-momentum-split.jld2")
    process_species(CR_e_gdf_momentum_filename, "electrons"; bins, outdir, should_title)

    return
end

function process_species(filename, species_name; bins, outdir, should_title = false)
    (log_p_nat, mle_fit_dists, hist_fit_dists) = dist_fits(filename, bins)
    @info("Fitted $species_name data to distributions")
    @info("Starting plots of $species_name")
    species_plots(
        log_p_nat, mle_fit_dists, hist_fit_dists;
        species_name, bins, outdir, should_title,
    )
    @info("Finished plots of $species_name")
end

function species_plots(
        log_p, mle_fit_dists, hist_fit_dists;
        species_name, bins, outdir, should_title,
    )
    fig = Figure()
    ax = Axis(fig[1, 1]; xticks = unit_ticks(log_p), axis_properties...)
    if should_title
        ax.title = plot_title
    end

    # cutoff line
    hlines!(ax, 2.0e-3, linewidth = 1.25, linestyle = :dash, color = :black)

    distmin, distmax = Inf, -Inf # for setting ticks
    valid_dist = filter(!isnan) ∘ skipmissing
    for (bin_count, dist) in zip(bins, hist_fit_dists)
        @debug("Processing bin count=$bin_count", dist)
        distances = bcdistances(mle_fit_dists.pf, dist.pf)
        distmin = minimum(distances |> valid_dist; init = distmin)
        distmax = maximum(distances |> valid_dist; init = distmax)
        linestyle = bin_count < 75 ? Linestyle([0,6,8]) : :solid # discriminator
        scatterlines!(ax, log_p, distances; label = "$bin_count", markersize, linestyle)
    end
    distmin = floor(Int, log10(distmin))
    distmax = ceil(Int, log10(distmax))
    ax.yticks = LogTicks(distmin:distmax)

    leg = Legend(
        fig[1,1], ax, "Number of bins";
        tellwidth = false, tellheight = false,
        orientation = :horizontal, valign = :top, nbanks = 2,
        margin = (10, 10, 10, 10),
        patchsize = (34, 20),   # make the lines bigger so the pattern shows
    )
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
    xlabel = L"$\log\, p$ ($m_\text{p} c$)",
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

        "--log-scale",
        Dict(:help => "Whether the y-axis should be in log-scale", :action => :store_true),

        "--title",
        Dict(:help => "Whether we should add a title to the figure", :action => :store_true),
    )
    return s
end
