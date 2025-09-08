using JLD2
using DataFrames
using Glob
using DelimitedFiles
using Format

using MCScatteringDataAnalysis


# A data row must not start with 3333..., those lines are for Fortran's pgf plotter
const data_row_predicate = !startswith("3333")

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

const datadirs = runpath .* [format("Seed-{:0>3}", s) for s in seeds]

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
        filename = joinpath(dir, filename)
        df = readdlm(filteredstream(filename; predicate), DataFrame; colspec)

        # tag each thing with the same tag, based on which directory they're in
        # so each df here has the same value in all rows, but it'll be different
        # for the next iteration
        insertcols!(df, 1, tagname=>tag)

        append!(bigdf, df)
    end

    return bigdf
end

function read_multiple_file_over_all_dirs(
    dirs::AbstractVector, filename_pattern, colspec::AbstractVector{<:ColumnSpecification};
    dirtags = seeds, filetags, dirtagname = :initial_seed, filetagname = :iteration,
    predicate = data_row_predicate)

    length(dirs) == length(dirtags) || throw(DimensionMismatch(
        "number of directories and directory tags must be equal"))

    bigdf = DataFrame()

    # iterate through all directories
    # (assume all these files are grouped by the "tags", i.e. each in a
    # different directory has a different tag)
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


function (@main)()
    # read and save the coupled spectra
    @info("Reading coupled spectra data")
    spectradf = read_one_file_over_all_dirs(
        datadirs, "mc_coupled_spectra.dat", coupled_spectra_cols)
    @info("Writing coupled spectra data")
    save_object("coupled-spectra.jld2", spectradf)

    # read and save the coupled weights
    @info("Reading coupled weights data")
    weightsdf = read_one_file_over_all_dirs(
        datadirs, "mc_coupled_wts.dat", coupled_weights_cols)
    @info("Writing coupled weights data")
    save_object("coupled-weights.jld2", weightsdf)

    # read and save the esc files
    @info("Reading dN/dp esc data")
    escdf = read_one_file_over_all_dirs(datadirs, "mc_dNdp_esc.dat", esc_cols)
    @info("Writing dN/dp esc data")
    save_object("dNdp-esc.jld2", escdf)

    # read and save the grid files
    @info("Reading grid data")
    griddf = read_one_file_over_all_dirs(datadirs, "mc_grid.dat", grid_cols)
    @info("Writing grid data")
    save_object("grid.jld2", griddf)


    # iterate through the dNdp on CR grid
    @info("Reading cosmic ray dN/dp data")
    CRdf = read_multiple_file_over_all_dirs(
        datadirs, "mc_dNdp_grid_CR_*.dat", CR_cols, filetags = 1:20)
    # save the dNdp on CR grid
    @info("Writing cosmic ray dN/dp data")
    save_object("dNdp-CR.jld2", CRdf)

    # iterate through the dNdp on therm grid
    @info("Reading thermal dN/dp data")
    thermdf = read_multiple_file_over_all_dirs(
        datadirs, "mc_dNdp_grid_therm_*.dat", therm_cols, filetags = 1:20)
    # save the dNdp on therm grid
    @info("Writing thermal dN/dp data")
    save_object("dNdp-therm.jld2", thermdf)
end
