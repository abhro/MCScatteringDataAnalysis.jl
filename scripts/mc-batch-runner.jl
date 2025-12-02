using DrWatson
@quickactivate "MCScatteringDataAnalysis"

using MCScatteringDataAnalysis.BatchProcessingUtilities: runmcproram, runpath, seeds
@main(args) = runmcprogram(string.(runpath, seeds))
