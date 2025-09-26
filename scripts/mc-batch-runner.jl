using DrWatson
@quickactivate "MCScatteringDataAnalysis"

using Dates

include("mc-batch-params.jl")

const timefmt = "e HH:MM:SS"
const progpath = `../../mc_cr.exe`

function runmcprogram(dirnames)
    for dirname in dirnames
        time = Dates.format(now(), timefmt)
        @info("Starting run in $dirname ($time)")

        cd(dirname) do
            # run mc.exe with mc_in piped
            pipeline(
                progpath, stdin = param_filename,
                stdout = "logs-stdout.txt", stderr = "logs-stderr.txt") |> run
        end
        time = Dates.format(now(), timefmt)
        @info("Ended run in $dirname ($time)")
    end
end

(@main)(args) = runmcprogram(string.(runpath, seeds))
