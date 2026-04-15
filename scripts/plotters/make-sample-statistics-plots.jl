using ArgParse: ArgParseSettings, add_arg_table!, parse_args
import Dates
using MCScatteringDataAnalysis: gdf_sample_stats, unit_ticks
using CairoMakie
using JLD2: load_object
using Statistics: mean, std
using StatsBase: StatsBase, skewness, kurtosis

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

    Makie.set_theme!(theme_latexfonts())
    Makie.update_theme!(
        fonts = Attributes(
            regular = "Libertinus Serif",
            bold = "Libertinus Serif Bold",
        )
    )

    CR_p_gdf_momentum_filename = joinpath(dir, "dNdp-CR-protons-momentum-split.jld2")
    CR_e_gdf_momentum_filename = joinpath(dir, "dNdp-CR-electrons-momentum-split.jld2")

    @info("Reading data ($(now()))")
    proton_stats = get_sample_stats(CR_p_gdf_momentum_filename; column)
    electron_stats = get_sample_stats(CR_e_gdf_momentum_filename; column)
    @info("Finished reading data ($(now()))")

    log_pₚ = proton_stats.p
    log_pₑ = electron_stats.p

    # for each sample statistic plot type

    # calculate sample statistics for protons and electrons
    p_means = proton_stats.means
    e_means = electron_stats.means
    p_std_devs = proton_stats.std_devs
    e_std_devs = electron_stats.std_devs

    # mean plot
    @info("Creating mean plots ($(now()))")
    make_sample_stat_plots(
        log_pₚ, log_pₑ, p_means, e_means;
        outdir, stat_title = "mean", ylabel = L"$⟨\log\,n_p⟩$"
    )

    @info("Creating std. dev. plots ($(now()))")
    make_std_dev_plots(log_pₚ, log_pₑ, p_std_devs, e_std_devs; outdir)

    # mean with std_dev envelope
    @info("Creating mean with std. dev. envelope plots ($(now()))")
    make_envelope_plots(log_pₚ, log_pₑ, p_means, e_means, p_std_devs, e_std_devs; outdir)

    # skewness
    @info("Creating skewness plots ($(now()))")
    make_sample_stat_plots(
        log_pₚ, log_pₑ, proton_stats.skewness, electron_stats.skewness;
        outdir, stat_title = "skewness", ylabel = L"γ"
    )

    # kurtosis
    @info("Creating kurtosis plots ($(now()))")
    make_sample_stat_plots(
        log_pₚ, log_pₑ, proton_stats.kurtosis, electron_stats.kurtosis;
        outdir, stat_title = "kurtosis", ylabel = "kurtosis"
    )

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

"""
Create proton, electron, and combined plots for a particular sample statistic

### Arguments
- `log_pₚ`: Vector containing (log of) proton momenta in mp*c units
- `log_pₑ`: Vector containing (log of) electron momenta in mp*c units
- `p_stat`: Sample statistic corresponding to each proton momentum
- `e_stat`: Sample statistic corresponding to each slice momentum
- `outdir`: Directory where the plots should be saved
- `stat_title`: Name of the sample statistic (mean/kurtosis/etc...)
- `ylabel`: Label for the plots' y-axes
"""
function make_sample_stat_plots(log_pₚ, log_pₑ, p_stat, e_stat; outdir, stat_title, ylabel)
    species_map = [
        ("protons", log_pₚ, p_stat, color_pf_p),
        ("electrons", log_pₑ, e_stat, color_pf_e),
    ]
    for (species, log_p, stat, color) in species_map
        @info("Making $species plot")
        fig, ax = make_figax(; stat_title, ylabel)
        lines!(ax, log_p, stat; color, label = species)
        ax.xticks = unit_ticks(log_p)
        #axislegend(ax) #, framevisible = false
        save(joinpath(outdir, "$species-$stat_title.svg"), fig)
    end

    # combined plot
    @info("Making combined plot")
    fig, ax = make_figax(; stat_title, ylabel)
    ax.xticks = unit_ticks(vcat(log_pₚ, log_pₑ))
    lines!(ax, log_pₚ, p_stat, color = color_pf_p, label = "protons")
    lines!(ax, log_pₑ, e_stat, color = color_pf_e, label = "electrons")
    axislegend(ax) #, framevisible = false
    save(joinpath(outdir, "combined-$stat_title.svg"), fig)
    return
end

function make_std_dev_plots(log_pₚ, log_pₑ, σₚ, σₑ; outdir)
    # plotting configs
    stat_title = "standard deviation"
    ylabel = L"σ"
    species_map = [
        ("protons", log_pₚ, σₚ, color_pf_p),
        ("electrons", log_pₑ, σₑ, color_pf_e),
    ]

    for (species, log_p, σ, color) in species_map
        @info("Making $species plot")
        fig, ax = make_figax(; stat_title, ylabel)
        ax.xticks = unit_ticks(log_p)
        lines!(ax, log_p, σ; color, label = species)
        #axislegend(ax, framevisible = false, position = :lt)
        save(joinpath(outdir, "$species-std-devs.svg"), fig)
        ax.yscale = log10
        save(joinpath(outdir, "$species-std-devs-logscale.svg"), fig)
    end

    # combined plot
    @info("Making combined plot")
    fig, ax = make_figax(; stat_title, ylabel)
    ax.xticks = unit_ticks(vcat(log_pₚ, log_pₑ))
    lines!(ax, log_pₚ, σₚ, color = color_pf_p, label = "protons")
    lines!(ax, log_pₑ, σₑ, color = color_pf_e, label = "electrons")
    axislegend(ax, position = :lt) #, framevisible = false
    save(joinpath(outdir, "combined-std-devs.svg"), fig)
    ax.yscale = log10
    save(joinpath(outdir, "combined-std-devs-logscale.svg"), fig)
    return
end

function make_envelope_plots(log_pₚ, log_pₑ, μₚ, μₑ, σₚ, σₑ; outdir)
    # plotting configs
    stat_title = "mean"
    ylabel = L"$⟨\log\,n_p⟩$"
    alpha = 0.4             # How transparent the error band should be
    species_map = [
        ("protons", log_pₚ, μₚ, σₚ, color_pf_p),
        ("electrons", log_pₑ, μₑ, σₑ, color_pf_e),
    ]
    for (species, log_p, μ, σ, color) in species_map
        @info("Making $species plot")
        fig, ax = make_figax(; stat_title, ylabel)
        ax.xticks = unit_ticks(log_p)
        lines!(ax, log_p, μ; color, label = species)
        band!(ax, log_p, μ + σ, μ - σ; alpha, color, label = species)
        #axislegend(ax, merge = true) #, framevisible = false
        save(joinpath(outdir, "$species-means-w-envelope.svg"), fig)
    end

    # combined plot
    @info("Making combined plot")
    fig, ax = make_figax(; stat_title, ylabel)
    ax.xticks = unit_ticks(vcat(log_pₚ, log_pₑ))
    lines!(ax, log_pₚ, μₚ, color = color_pf_p, label = "protons")
    band!(ax, log_pₚ, μₚ + σₚ, μₚ - σₚ; alpha, color = color_pf_p, label = "protons")
    lines!(ax, log_pₑ, μₑ, color = color_pf_e, label = "electrons")
    band!(ax, log_pₑ, μₑ + σₑ, μₑ - σₑ; alpha, color = color_pf_e, label = "electrons")
    axislegend(ax, merge = true) #, framevisible = false
    save(joinpath(outdir, "combined-means-w-envelope.svg"), fig)
    return
end

const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors();
const axis_properties = (;
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorticksvisible = true,
    yminorticksvisible = true,
    xlabel = L"$\log\,p$ $(m_\text{p} c)$",
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
