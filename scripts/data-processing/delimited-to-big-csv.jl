using DrWatson
@quickactivate "MCScatteringDataAnalysis"

using CSV
using DataFrames: AbstractDataFrame, DataFrame, insertcols!
using Glob: glob
using DelimitedFiles
using Format: format
import Comonicon

using MCScatteringDataAnalysis
using MCScatteringDataAnalysis: ColumnSpecification, name, uses_sentinels,
                                coupled_spectra_cols, coupled_weights_cols,
                                esc_cols, grid_cols, CR_cols, therm_cols

function (@main)(args = [])
    if length(args) != 2
        println("Usage: $PROGRAM_FILE <raw data directory> <processed data directory>")
        println("Note: processed data directory must be relative to the project data path")
        exit(1)
    end
    datadirs = glob("Seed-*", args[1])
    tags = parse.(Int, replace.(datadirs, r".*Seed-0{0,3}" => ""))

    outdir = datadir(args[2])

    # read and save the coupled spectra
    @info("Reading coupled spectra data")
    spectradf = read_one_file_over_all_dirs(
        datadirs, "mc_coupled_spectra.dat", coupled_spectra_cols; tags
    )
    outfilepath = joinpath(outdir, "coupled-spectra.csv.gz")
    @info("Writing coupled spectra data to $outfilepath")
    CSV.write(outfilepath, spectradf; compress = true)

    # read and save the coupled weights
    @info("Reading coupled weights data")
    weightsdf = read_one_file_over_all_dirs(
        datadirs, "mc_coupled_wts.dat", coupled_weights_cols; tags
    )
    outfilepath = joinpath(outdir, "coupled-weights.csv.gz")
    @info("Writing coupled weights data to $outfilepath")
    CSV.write(outfilepath, weightsdf; compress = true)

    # read and save the esc files
    @info("Reading dN/dp esc data")
    escdf = read_one_file_over_all_dirs(datadirs, "mc_dNdp_esc.dat", esc_cols; tags)
    outfilepath = joinpath(outdir, "dNdp-esc.csv.gz")
    @info("Writing dN/dp esc data to $outfilepath")
    CSV.write(outfilepath, escdf; compress = true)

    # read and save the grid files
    @info("Reading grid data")
    griddf = read_one_file_over_all_dirs(datadirs, "mc_grid.dat", grid_cols; tags)
    @info("Writing grid data to $outfilepath")
    outfilepath = joinpath(outdir, "grid.csv.gz")
    CSV.write(outfilepath, griddf; compress = true)

    # iterate through the dNdp on CR grid
    @info("Reading cosmic ray dN/dp data")
    CRdf = read_multiple_file_over_all_dirs(
        datadirs, "mc_dNdp_grid_CR_*.dat", CR_cols; filetags = 1:20, dirtags = tags
    )
    # save the dNdp on CR grid
    outfilepath = joinpath(outdir, "dNdp-CR.csv.gz")
    @info("Writing cosmic ray dN/dp data to $outfilepath")
    CSV.write(outfilepath, CRdf; compress = true)

    # iterate through the dNdp on therm grid
    @info("Reading thermal dN/dp data")
    thermdf = read_multiple_file_over_all_dirs(
        datadirs, "mc_dNdp_grid_therm_*.dat", therm_cols; filetags = 1:20, dirtags = tags
    )
    # save the dNdp on therm grid
    outfilepath = joinpath(outdir, "dNdp-therm.csv.gz")
    @info("Writing thermal dN/dp data to $outfilepath")
    CSV.write(outfilepath, thermdf; compress = true)

    return 0
end
