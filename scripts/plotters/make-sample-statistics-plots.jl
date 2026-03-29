using ArgParse: ArgParseSettings, add_arg_table!, parse_args
import Dates
using MCScatteringDataAnalysis: gdf_sample_stats
using CairoMakie
using JLD2: load_object
using Statistics: mean, std
using StatsBase: StatsBase, skewness, kurtosis

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

    # get option for if we want combined, protons only, and/or electrons only
    #plot_combined = args["combined"]
    #plot_protons = args["protons"]
    #plot_electrons = args["electrons"]

    #println(args)
    #if !plot_combined && !plot_protons && !plot_electrons
    #    println(stderr, "Please select at least one configuration to plot. (Use -h/--help for help)")
    #    exit(1)
    #end

    CR_p_gdf_momentum_filename = joinpath(dir, "dNdp-CR-protons-momentum-split.jld2")
    CR_e_gdf_momentum_filename = joinpath(dir, "dNdp-CR-electrons-momentum-split.jld2")

    @info("Reading data ($(now()))")
    proton_stats = get_sample_stats(CR_p_gdf_momentum_filename; column)
    electron_stats = get_sample_stats(CR_e_gdf_momentum_filename; column)
    @info("Finished reading data ($(now()))")

    proton_log_p_nat = proton_stats.p
    electron_log_p_nat = electron_stats.p

    # for each sample statistic plot type

    # calculate sample statistics for protons and electrons
    p_means = proton_stats.means
    e_means = electron_stats.means
    p_std_devs = proton_stats.std_devs
    e_std_devs = electron_stats.std_devs

    # mean plot
    @info("Creating mean plots ($(now()))")
    # proton plot
    @info("Making proton plot")
    fig, ax = make_figax(stat_title = "mean", ylabel = "⟨log nₚ⟩")
    lines!(ax, proton_log_p_nat, p_means, color = color_pf_p, label = "protons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "proton-means.svg"), fig)

    # electron plot
    @info("Making electron plot")
    fig, ax = make_figax(stat_title = "mean", ylabel = "⟨log nₚ⟩")
    lines!(ax, electron_log_p_nat, e_means, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "electron-means.svg"), fig)

    # combined plot
    @info("Making combined plot")
    fig, ax = make_figax(stat_title = "mean", ylabel = "⟨log nₚ⟩")
    lines!(ax, proton_log_p_nat, p_means, color = color_pf_p, label = "protons")
    lines!(ax, electron_log_p_nat, e_means, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "combined-means.svg"), fig)

    @info("Creating std_dev plots ($(now()))")
    # proton plot
    @info("Making proton plot")
    fig, ax = make_figax(stat_title = "standard deviation", ylabel = "σ")
    lines!(ax, proton_log_p_nat, p_std_devs, color = color_pf_p, label = "protons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "proton-std-devs.svg"), fig)
    ax.yscale = log10
    save(joinpath(outdir, "proton-std-devs-logscale.svg"), fig)

    # electron plot
    @info("Making electron plot")
    fig, ax = make_figax(stat_title = "standard deviation", ylabel = "σ")
    lines!(ax, electron_log_p_nat, e_std_devs, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "electron-std-devs.svg"), fig)
    ax.yscale = log10
    save(joinpath(outdir, "electron-std-devs-logscale.svg"), fig)

    # combined plot
    @info("Making combined plot")
    fig, ax = make_figax(stat_title = "standard deviation", ylabel = "σ")
    lines!(ax, proton_log_p_nat, p_std_devs, color = color_pf_p, label = "protons")
    lines!(ax, electron_log_p_nat, e_std_devs, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "combined-std-devs.svg"), fig)
    ax.yscale = log10
    save(joinpath(outdir, "combined-std-devs-logscale.svg"), fig)


    # mean with std_dev envelope
    @info("Creating mean with std_dev envelope plots ($(now()))")
    # proton plot
    @info("Making proton plot")
    fig, ax = make_figax(stat_title = "mean", ylabel = "⟨log nₚ⟩")
    lines!(ax, proton_log_p_nat, p_means, color = color_pf_p, label = "protons")
    band!(ax, proton_log_p_nat, p_means + p_std_devs, p_means - p_std_devs, alpha = 0.4, color = color_pf_p, label = "protons")
    axislegend(ax, framevisible = false, merge = true)
    save(joinpath(outdir, "proton-means-w-envelope.svg"), fig)

    # electron plot
    @info("Making electron plot")
    fig, ax = make_figax(stat_title = "mean", ylabel = "⟨log nₚ⟩")
    lines!(ax, electron_log_p_nat, e_means, color = color_pf_e, label = "electrons")
    band!(ax, electron_log_p_nat, e_means + e_std_devs, e_means - e_std_devs, alpha = 0.4, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false, merge = true)
    save(joinpath(outdir, "electron-means-w-envelope.svg"), fig)

    # combined plot
    @info("Making combined plot")
    fig, ax = make_figax(stat_title = "mean", ylabel = "⟨log nₚ⟩")
    lines!(ax, proton_log_p_nat, p_means, color = color_pf_p, label = "protons")
    band!(ax, proton_log_p_nat, p_means + p_std_devs, p_means - p_std_devs, alpha = 0.4, color = color_pf_p, label = "protons")
    lines!(ax, electron_log_p_nat, e_means, color = color_pf_e, label = "electrons")
    band!(ax, electron_log_p_nat, e_means + e_std_devs, e_means - e_std_devs, alpha = 0.4, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false, merge = true)
    save(joinpath(outdir, "combined-means-w-envelope.svg"), fig)

    # skewness
    @info("Creating skewness plots ($(now()))")
    make_skewness_plots(proton_log_p_nat, electron_log_p_nat, proton_stats.skewness, electron_stats.skewness; column, outdir)

    # kurtosis
    @info("Creating kurtosis plots ($(now()))")
    make_kurtosis_plots(proton_log_p_nat, electron_log_p_nat, proton_stats.kurtosis, electron_stats.kurtosis; column, outdir)

    return
end

now() = Dates.format(Dates.now(), "HH:MM:SS")

function get_sample_stats(filename; column)
    gdf = load_object(filename)
    p = keys(gdf) .|> values .|> first
    means = gdf_sample_stats(mean, gdf; column)
    std_devs = gdf_sample_stats(std, gdf; column)
    skewness = gdf_sample_stats(StatsBase.skewness, gdf; column)
    kurtosis = gdf_sample_stats(StatsBase.kurtosis, gdf; column)

    return (; p, means, std_devs, skewness, kurtosis)
end

function make_figax(; stat_title, ylabel)
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Sample $stat_title vs momentum slice",
        axis_properties...,
        ylabel,
    )
    return fig, ax
end

function make_skewness_plots(proton_log_p_nat, electron_log_p_nat, p_skewness, e_skewness; column, outdir)
    # proton plot
    @info("Making proton plot")
    fig, ax = make_figax(stat_title = "skewness", ylabel = "γ")
    lines!(ax, proton_log_p_nat, p_skewness, color = color_pf_p, label = "protons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "proton-skewness.svg"), fig)

    # electron plot
    @info("Making electron plot")
    fig, ax = make_figax(stat_title = "skewness", ylabel = "γ")
    lines!(ax, electron_log_p_nat, e_skewness, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "electron-skewness.svg"), fig)

    # combined plot
    @info("Making combined plot")
    fig, ax = make_figax(stat_title = "skewness", ylabel = "γ")
    lines!(ax, proton_log_p_nat, p_skewness, color = color_pf_p, label = "protons")
    lines!(ax, electron_log_p_nat, e_skewness, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "combined-skewness.svg"), fig)

    return
end

function make_kurtosis_plots(proton_log_p_nat, electron_log_p_nat, p_kurtosis, e_kurtosis; column, outdir)
    # proton plot
    @info("Making proton plot")
    fig, ax = make_figax(stat_title = "kurtosis", ylabel = "kurtosis")
    lines!(ax, proton_log_p_nat, p_kurtosis, color = color_pf_p, label = "protons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "proton-kurtosis.svg"), fig)

    # electron plot
    @info("Making electron plot")
    fig, ax = make_figax(stat_title = "kurtosis", ylabel = "kurtosis")
    lines!(ax, electron_log_p_nat, e_kurtosis, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "electron-kurtosis.svg"), fig)

    # combined plot
    @info("Making combined plot")
    fig, ax = make_figax(stat_title = "kurtosis", ylabel = "kurtosis")
    lines!(ax, proton_log_p_nat, p_kurtosis, color = color_pf_p, label = "protons")
    lines!(ax, electron_log_p_nat, e_kurtosis, color = color_pf_e, label = "electrons")
    axislegend(ax, framevisible = false)
    save(joinpath(outdir, "combined-kurtosis.svg"), fig)

    return
end

const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();
const axis_properties = (;
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorticksvisible = true,
    yminorticksvisible = true,
    xlabel = "log p (nat)",
)

"""
Create and return an ArgParse.jl arg table
"""
function get_parser()
    s = ArgParseSettings(; help_width = 108)

    add_arg_table!(
        s,
        "data-dir",
        Dict(:help => "Where to look for the .jld2 files", :required => true),
        "out-dir",
        Dict(:help => "Where to save the plots", :required => true),

        "--column",
        Dict(:help => "Which column of the dataframe to take statistics off", :default => :log_dNdp_cr_pf, :metavar => "COL"),

        #"--no-combined", Dict(:dest_name => "combined", :action => :store_false),
        #"--combined",
        #Dict(:help => "Whether we should plot both proton and electron statistics on the same graph", :action => :store_true),

        #"--no-protons", Dict(:dest_name => "protons", :action => :store_false),
        #"--protons",
        #Dict(:help => "Whether we should plot protons on its own graph", :action => :store_true),

        #"--no-electrons", Dict(:dest_name => "electrons", :action => :store_false),
        #"--electrons",
        #Dict(:help => "Whether we should plot electrons on its own graph", :action => :store_true),
    )
    return s
end
