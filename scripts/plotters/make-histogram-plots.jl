using JLD2: load_object
using DataFrames
using LinearAlgebra: normalize
using CairoMakie
using StatsBase
using Distributions: Normal, params
using Printf: @sprintf
using ArgParse: ArgParseSettings, add_arg_table!, parse_args
using MCScatteringDataAnalysis: fit_dist_to_histogram

const normalization = :pdf

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
    should_plot_mle_fit = args["mle-fit"]
    should_plot_hist_fit = args["hist-fit"]

    CR_p_gdf_momentum_filename = joinpath(dir, "dNdp-CR-protons-momentum-split.jld2")
    CR_p_gdf_momentum = load_object(CR_p_gdf_momentum_filename)
    proton_log_p_nat = keys(CR_p_gdf_momentum) .|> values .|> first

    @info("Starting plots of protons")
    process_species(CR_p_gdf_momentum, proton_log_p_nat, "proton", outdir;
                    column, bins, yscale, should_plot_mle_fit, should_plot_hist_fit)
    @info("Finished plots of protons")

    if args["plot-electrons"]
        CR_e_gdf_momentum_filename = joinpath(dir, "dNdp-CR-electrons-momentum-split.jld2")
        CR_e_gdf_momentum = load_object(CR_e_gdf_momentum_filename)
        electron_log_p_nat = keys(CR_e_gdf_momentum) .|> values .|> first

        @info("Starting plots of electrons")
        process_species(CR_e_gdf_momentum, electron_log_p_nat, "electron", outdir;
                        column, bins, yscale, should_plot_mle_fit, should_plot_hist_fit)
        @info("Finished plots of electrons")
    end

    return
end

function process_species(
        species_gdf, log_p_all, species_name, outdir;
        column, bins, yscale, should_plot_mle_fit, should_plot_hist_fit,
    )
    # loop over momentum slices in particle species
    for (log_p, df) in zip(log_p_all, species_gdf)
        log_nₚ = df[!, column] |> skipmissing |> collect

        if isempty(log_nₚ)
            continue
        end

        # make histogram of the thing -- ends up not being used :/
        ##histfit = fit(Histogram, log_nₚ; nbins = bins)
        ##histfit = normalize(histfit, mode = :pdf)
        histfit = fit_dist_to_histogram(Normal, log_nₚ; nbins = bins)

        # make normal distribution of the thing
        distfit = fit(Normal, log_nₚ)

        # make the plot of the thing
        fig = make_plot(log_nₚ, distfit, histfit; bins, yscale, should_plot_mle_fit, should_plot_hist_fit)
        filename = @sprintf("%s-momentum-slice-histogram-log-p-%.01f-mpc.svg", species_name, log_p)
        save(joinpath(outdir, filename), fig)
    end
end

const axis_properties = (;
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorticksvisible = true,
    yminorticksvisible = true,
    xlabel = "log(nₚ)",
)

const color_pf_p, color_sf_p, color_ISM_p, color_pf_e, color_sf_e, color_ISM_e = Makie.wong_colors()

"""
Plot a histogram of samples and associated distribution.
"""
function make_plot(samples, mle_dist, hist_dist; bins, yscale, should_plot_mle_fit, should_plot_hist_fit)
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        ylabel = "pdf",
        yscale,
        ##title = "Histogram of protons dN/dp at log p = $log_p_nat_at_slice_p (mₚc)",
        axis_properties...
    )

    n = length(samples)
    ##samples ./= std(samples)
    if n != 0
        stephist!(ax, samples; label = "plasma frame ($n samples)", bins, normalization, color = color_pf_p, linewidth = 2)
    end

    if should_plot_mle_fit && !ismissing(mle_dist)
        # n points seems like a good heuristic
        data_span = range(extrema(samples)..., length = max(n, 2000))
        dist_label = @sprintf("MLE fit (μ=%.3f, σ=%.3f)", params(mle_dist)...)
        plot!(ax, data_span, mle_dist, label = dist_label, color = :indianred, linewidth = 1)
    end
    if should_plot_hist_fit && !ismissing(hist_dist)
        # n points seems like a good heuristic
        data_span = range(extrema(samples)..., length = max(n, 2000))
        dist_label = @sprintf("Curve fit (μ=%.3f, σ=%.3f)", params(hist_dist)...)
        plot!(ax, data_span, hist_dist, label = dist_label, color = :brown, linestyle = :dash)
    end

    axislegend(ax, framevisible = false)

    return fig
end

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

        "--plot-electrons",
        Dict(:help => "Whether we should also plot electrons", :action => :store_true),
        "--no-plot-electrons",
        Dict(:dest_name => "plot-electrons", :action => :store_false),

        "--mle-fit",
        Dict(
            :help => "Whether to overlay a plot of the MLE fitted Gaussian distribution",
            :action => :store_true,
        ),
        "--no-mle-fit",
        Dict(:dest_name => "mle-fit", :action => :store_false),

        "--hist-fit",
        Dict(
            :help => "Whether to overlay a plot of the histogram fitted Gaussian distribution",
            :action => :store_true,
        ),
        "--no-hist-fit",
        Dict(:dest_name => "hist-fit", :action => :store_false),

        "--log-scale",
        Dict(:help => "Whether the y-axis should be in log-scale", :action => :store_true),
    )
    return s
end
