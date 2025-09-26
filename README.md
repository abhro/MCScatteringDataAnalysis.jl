# MCScatteringDataAnalysis.jl

This code base is using the [Julia Language](https://julialang.org/) and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project.

This project contains of scripts and notebooks for analyzing the data created by the MonteCarlo\_cr program.

For the original Fortran program producing the data, see the [MonteCarlo\_cr repository](https://github.fit.edu/dwarren/MonteCarlo_cr).

## Setting up the environment
Ensure you have [Julia installed](https://julialang.org/install) on your device. Then, clone the repository.
Notice that raw data are typically not included in the git-history and may need to be downloaded independently. Please contact the repository authors/maintainers for a copy of the raw data.

Navigating to the root directory of the cloned repository, open a Julia REPL. Within the Julia REPL, run
```julia
julia> ]            # Press `]` to drop into pkg-mode
pkg> add DrWatson   # Install DrWatson.jl globally, for using `quickactivate`
   Resolving package versions...
[...]

pkg> activate .
  Activating project at [...]

(MCScatteringDataAnalysis) pkg> instantiate
[...]
```

Alternatively, you can run directly from the command line
```sh
$ julia -e 'import Pkg; Pkg.add("DrWatson")'
$ julia --project -e 'import Pkg; Pkg.instantiate()'
```

This will install all the necessary (prerequisite) packages for you to be able to run the scripts and everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "MCScatteringDataAnalysis"
```
which auto-activate the project and enable local path handling from DrWatson.

## Running the notebooks
Make sure the project environment has been instantiated. Then, from the project root directory, run
```sh
$ julia --project notebooks/runpluto.jl
```
If you want to also open up the browser, append `--browser` to the end of the command line string.
