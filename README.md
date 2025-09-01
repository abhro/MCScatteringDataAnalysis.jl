# MCScatteringDataAnalysis.jl

Collection of scripts and notebooks for analyzing the data created by the MonteCarlo\_cr program.

For the original Fortran program producing the data, see the [MonteCarlo\_cr repository](https://github.fit.edu/dwarren/MonteCarlo_cr).

## Setting up the environment
Ensure you have [Julia installed](https://julialang.org/install) on your device. Then, clone the repository.

Navigating to the root directory of the cloned repository, open a Julia REPL. Within the Julia REPL, run
```julia
julia> ] # Press `]` to drop into pkg mode
pkg> activate .
  Activating project at [...]

(MCScatteringDataAnalysis) pkg> instantiate
[...]

(MCScatteringDataAnalysis) pkg> activate notebooks
  Activating project at [...]

(notebooks) pkg> instantiate
[...]
```

This will install all the necessary prerequisites. Alternatively, you can run directly from the command line
```sh
$ julia --project -e 'import Pkg; Pkg.instantiate()'
$ julia --project=notebooks -e 'import Pkg; Pkg.instantiate()'
```

## Running the notebooks
For the notebooks, make sure the project environment and the enclosed `notebooks` environment has been instantiated. Then in the `notebooks` directory, run
```sh
$ julia --project runpluto.jl
```
Alternatively, from the project root directory, run
```sh
$ julia --project=notebooks notebooks/runpluto.jl
```
If you want to also open up the browser, append `--browser` to the end of the command line string.
