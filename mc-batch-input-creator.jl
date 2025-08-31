#!julia --project

using Format

include("mc-batch-params.jl")

"""
    writeparams(filename, seed)

Write Fortran parameter block using `params` from mc-batch-params.jl to `filename`.
Includes `seed` as the `ISEED` parameter
"""
function writeparams(filename, seed)
    seedline = "ISEED $seed Seed for random number generator"
    open(filename, "w") do f
        println(f, seedline)
        println(f, params)
    end
end

function initdirs(seeds)
    # TODO use 3 digit seed names
    dirnames = string.(runpathbase, seeds)
    for (seed, dirname) in zip(seeds, dirnames)
        # Create the folder
        mkdir(dirname)

        # cd into the folder
        cd(dirname) do
            writeparams(param_filename, seed)
        end
        @info("Created parameter file for seed = $seed at $dirname")
    end
    return dirnames
end

initdirs(seeds)
