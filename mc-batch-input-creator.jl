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
    mkpath(runpath)
    cd(runpath)
    for seed in seeds
        dirname = format("Seed-{:0>3}", seed)
        # Create the folder
        mkdir(dirname)

        # `cd()` into the folder and write the parameters
        cd(dirname) do
            writeparams(param_filename, seed)
            @info("Created parameter file for seed = $seed at $(pwd())")
        end
    end
end

initdirs(seeds)
