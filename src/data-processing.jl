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
A DataFrame containing.... (TODO: finish writing docstring)
"""
function read_one_file_over_all_dirs(
        dirs::AbstractVector, filename, colspec::AbstractVector{<:ColumnSpecification};
        tags = seeds, tagname = :initial_seed, predicate = data_row_predicate
    )

    if length(dirs) != length(tags)
        throw(DimensionMismatch("number of directories and tags must be equal"))
    end

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
        insertcols!(df, 1, tagname => tag)

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
        predicate = data_row_predicate
    )

    if length(dirs) != length(dirtags)
        throw(DimensionMismatch("number of directories and directory tags must be equal"))
    end

    bigdf = DataFrame()

    # iterate through all directories
    # (assume all these files are grouped by the "tags",
    # i.e., each in a different directory has a different tag)
    for (dirtag, dir) in zip(dirtags, dirs)

        # go in to the directory first before reading multiple files
        cd(dir) do
            filenames = glob(filename_pattern) |> sort
            if length(filenames) != length(filetags)
                throw(
                    DimensionMismatch(
                        "In directory $dir: number of files must match number of file tags"
                    )
                )
            end

            for (filetag, filename) in zip(filetags, filenames)
                df = readdlm(filteredstream(filename; predicate), DataFrame; colspec)

                insertcols!(df, 1, dirtagname => dirtag, filetagname => filetag)

                append!(bigdf, df)
            end
        end
    end

    return bigdf
end


"""
De-histogram each `DataFrame` of a `GroupedDataFrame`, then return one big DataFrame
containing all the rows.
"""
function dehistogram(gdf::GroupedDataFrame)
    df_dehistogrammed = DataFrame()
    for df in gdf
        # de-histogram each of the sub-dataframes (make a copy first)
        df_dehist = df[(begin + 1):2:end, :]

        # XXX XXX XXX XXX
        # **Need to explore `psd_mom_bounds`**, and how data repeats in adjacent cells (high/low bounds)
        # -> Choose the upper bounds for de-histogramming the data, because the lowest
        #    bin covers -99.0 to -19.3 g⋅cm/s, which is too wide. The highest bin covers
        #    -2.3997 to -2.2997 g⋅cm/s
        # -> runs into a fencepost problem. complete last cell (343 for each grouped df)
        #    has a unrepeated value?

        # push that copy into a bigger temp df
        # the cols attribute is an additional sanity check
        append!(df_dehistogrammed, df_dehist, cols = :orderequal)
    end
    return df_dehistogrammed
end
