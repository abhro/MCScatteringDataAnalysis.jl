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

## Processing the raw data
For each distinct mc\_in.txt file that's being fed to the MonteCarlo\_cr Fortran program, put each of them in a distinct folder. You may group these folders as is convenient. Generally, you do not want to mix folders that contain different physical parameters, but only differ in the initial RNG seed. Here, the common approach is to create a folder for each set of distinct physical parameters, then within that folder, create subfolders titled as Seed-NNNN where NNNN is a 4-digit number from 0000 to 9999. Once you're in the top-level directory, you can aggregate all the different simulation data using the `delimited-to-hdf5.jl` script. (This actually no longer creates HDF5 files, so the script should be renamed.) For example,
```sh
$ ls data/physical-parameters-whatever/.../
Seed-0000   Seeed-0001   Seed-0002    [...]
$ julia --project scripts/delimited-to-hdf5.jl data/physical-parameters-whatever/.../ physical-parameters-whatever-processed
```

This will create a [gzipped](https://gzip.org/) CSV file for each of the datasets the Fortran program outputs. Namely, the coupled spectra and weights data, the dN/dp data for thermal, CR (cosmic ray) and escaping particles, and the grid data. After running the `delimited-to-hdf5.jl` script, the processed data folder should look something like this
```sh
$ ls data/ physical-parameters-whatever-processed
coupled-spectra.csv.gz  dNdp-CR-electrons.csv.gz  dNdp-CR.csv.gz   dNdp-therm.csv.gz
coupled-weights.csv.gz  dNdp-CR-protons.csv.gz    dNdp-esc.csv.gz  grid.csv.gz
```

## Running the notebooks
Make sure the project environment has been instantiated. Then, from the project root directory, run
```sh
$ julia --project notebooks/runpluto.jl
```
If you want to also open up the browser, append `--browser` to the end of the command line string.
