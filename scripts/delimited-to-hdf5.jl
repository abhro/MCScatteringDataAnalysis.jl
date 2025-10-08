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



# A data row must not start with 3333..., those lines are for Fortran's pgf plotter
const data_row_predicate = !startswith("3333")

"""
    filteredstream(filename::AbstractString; predicate = data_row_predicate)

Read `filename` and return a stream over it such that each line satisfies `predicate`.

### Arguments
- `filename`: Path to a file
- `predicate`: (keyword, optional) A function taking a `String` and
  returning a `Bool` that determines the filter to be applied over each line.

### Returns
- `filtered_buffer`: An `IOBuffer` containing lines that match `predicate`
  (i.e., filters lines that fail `predicate`.)
"""
function filteredstream(filename::AbstractString; predicate = data_row_predicate)
    filtered_buffer = IOBuffer()
    infilestream = open(filename)
    for line in eachline(infilestream)
        if predicate(line)
            write(filtered_buffer, line)
            write(filtered_buffer, "\n")
        end
    end
    close(infilestream)

    seekstart(filtered_buffer)
    return filtered_buffer
end

# include file to get `seeds` and `runpath`
include("mc-batch-params.jl")

const datadirs = joinpath.(runpath, format.("Seed-{:0>4}", seeds))

"""
    DelimitedFiles.readdlm(source, T; colspec)

Type-pirated version of `readdlm()` that reads a delimited file into a DataFrame.
Note that `T` must be a type such that `T <: DataFrames.AbstractDataFrame`.

The specification for each column is given through `colspec`.
"""
function DelimitedFiles.readdlm(source, T::Type{<:AbstractDataFrame}; colspec)
    rawmat = readdlm(source)
    colnames = name.(colspec)
    df = T(rawmat, colnames)
    for (colname, column) in zip(colnames, colspec)
        newcol = convert.(eltype(column), df[!, colname])
        if uses_sentinels(column)
            # TODO allow custom sentinels
            # not using replace! because it changes the array type
            newcol = replace(newcol, -99 => missing)
        end
        df[!, colname] = newcol
    end

    return df
end

"""
    read_one_file_over_all_dirs(dirs, filename, colspec; tags, tagname, predicate)

### Arguments
- `dirs`
- `filename`
- `colspec`
- `tags`
- `tagname`
- `predicate`

### Returns
"""
function read_one_file_over_all_dirs(
    dirs::AbstractVector, filename, colspec::AbstractVector{<:ColumnSpecification};
    tags = seeds, tagname = :initial_seed, predicate = data_row_predicate)

    length(dirs) == length(tags) || throw(DimensionMismatch(
        "number of directories and tags must be equal"))

    bigdf = DataFrame()

    # iterate through all directories
    # (assume all these files are grouped by the "tags", i.e. each in a
    # different directory has a different tag)
    for (tag, dir) in zip(tags, dirs)
        filepath = joinpath(dir, filename)
        df = readdlm(filteredstream(filepath; predicate), DataFrame; colspec)

        # tag each thing with the same tag, based on which directory they're in
        # so each df here has the same value in all rows, but it'll be different
        # for the next iteration
        insertcols!(df, 1, tagname=>tag)

        append!(bigdf, df)
    end

    return bigdf
end

"""
    read_multiple_file_over_all_dirs(
        dirs, filename_pattern, colspec;
        dirtags, filetags, dirtagname, filetagname, predicate)

### Arguments
- `dirs`: List of directories to read files from.
- `filename_pattern`: Glob pattern for filenames within each directory
- `colspec`: Column specification to use when reading data file.
- `dirtags`: How to tag files using directory-level information.
- `filetags`: How to tag each file from the same directory.
- `dirtagname`:
- `filetagname`:
- `predicate`:

### Returns

A `DataFrame` containing columns as specified in `colspec`.
It also has two extra columns as specified by `dirtagname` and `filetagname`.
"""
function read_multiple_file_over_all_dirs(
    dirs::AbstractVector, filename_pattern, colspec::AbstractVector{<:ColumnSpecification};
    dirtags = seeds, filetags, dirtagname = :initial_seed, filetagname = :iteration,
    predicate = data_row_predicate)

    length(dirs) == length(dirtags) || throw(DimensionMismatch(
        "number of directories and directory tags must be equal"))

    bigdf = DataFrame()

    # iterate through all directories
    # (assume all these files are grouped by the "tags",
    # i.e., each in a different directory has a different tag)
    for (dirtag, dir) in zip(dirtags, dirs)

        # go in to the directory first before reading multiple files
        cd(dir) do
            filenames = glob(filename_pattern) |> sort
            length(filenames) == length(filetags) || throw(DimensionMismatch(
                "In directory $dir: number of files must match number of file tags"
            ))

            for (filetag, filename) in zip(filetags, filenames)
                df = readdlm(filteredstream(filename; predicate), DataFrame; colspec)

                insertcols!(df, 1, dirtagname => dirtag, filetagname => filetag)

                append!(bigdf, df)
            end
        end
    end

    return bigdf
end

function (@main)(args)
    datadirs = glob("Seed-*", runpath)
    replaced = replace.(datadirs, r".*Seed-0{0,3}"=>"")
    tags = parse.(Int, replaced)
    dirtags = tags

    outdir = args[1]

    # read and save the coupled spectra
    @info("Reading coupled spectra data")
    spectradf = read_one_file_over_all_dirs(
        datadirs, "mc_coupled_spectra.dat", coupled_spectra_cols; tags)
    outfilepath = joinpath(outdir, "coupled-spectra.csv.gz")
    @info("Writing coupled spectra data to $outfilepath")
    CSV.write(outfilepath, spectradf; compress=true)

    # read and save the coupled weights
    @info("Reading coupled weights data")
    weightsdf = read_one_file_over_all_dirs(
        datadirs, "mc_coupled_wts.dat", coupled_weights_cols; tags)
    outfilepath = joinpath(outdir, "coupled-weights.csv.gz")
    @info("Writing coupled weights data to $outfilepath")
    CSV.write(outfilepath, weightsdf; compress=true)

    # read and save the esc files
    @info("Reading dN/dp esc data")
    escdf = read_one_file_over_all_dirs(datadirs, "mc_dNdp_esc.dat", esc_cols; tags)
    outfilepath = joinpath(outdir, "dNdp-esc.csv.gz")
    @info("Writing dN/dp esc data to $outfilepath")
    CSV.write(outfilepath, escdf; compress=true)

    # read and save the grid files
    @info("Reading grid data")
    griddf = read_one_file_over_all_dirs(datadirs, "mc_grid.dat", grid_cols; tags)
    @info("Writing grid data to $outfilepath")
    outfilepath = joinpath(outdir, "grid.csv.gz")
    CSV.write(outfilepath, griddf; compress=true)

    # iterate through the dNdp on CR grid
    @info("Reading cosmic ray dN/dp data")
    CRdf = read_multiple_file_over_all_dirs(
        datadirs, "mc_dNdp_grid_CR_*.dat", CR_cols; filetags = 1:20, dirtags = tags)
    # save the dNdp on CR grid
    outfilepath = joinpath(outdir, "dNdp-CR.csv.gz")
    @info("Writing cosmic ray dN/dp data to $outfilepath")
    CSV.write(outfilepath, CRdf; compress=true)

    # iterate through the dNdp on therm grid
    @info("Reading thermal dN/dp data")
    thermdf = read_multiple_file_over_all_dirs(
        datadirs, "mc_dNdp_grid_therm_*.dat", therm_cols; filetags = 1:20, dirtags = tags)
    # save the dNdp on therm grid
    outfilepath = joinpath(outdir, "dNdp-therm.csv.gz")
    @info("Writing thermal dN/dp data to $outfilepath")
    CSV.write(outfilepath, thermdf; compress=true)
end
