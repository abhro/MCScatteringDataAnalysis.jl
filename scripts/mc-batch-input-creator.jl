using DrWatson
@quickactivate "MCScatteringDataAnalysis"

using MCScatteringDataAnalysis.BatchProcessingUtilities: initdirs
@main(args = []) = initdirs(seeds)
