# MCScatteringDataAnalysis.jl

Collection of scripts and notebooks for analyzing the data created by the MonteCarlo\_cr program.

## Setting up the environment
First, clone the repository. Then, in a Julia REPL, run
```julia
julia> ] # Press `]` to drop into pkg mode
pkg> activate .

(MCScatteringDataAnalysis) pkg> instantiate
```

This will install all the necessary prerequisites. Alternatively, you can run directly from the command line
```sh
$ julia --project -e 'import Pkg; Pkg.instantiate()'
```

## Running the notebooks
For the notebooks, make sure the project environment has been instantiated. Then in the project root, run
```sh
$ julia --project runpluto.jl
```
If you want to also open up the browser, append `--browser` to the end of the command line string.
