module BatchProcessingUtilities
using DrWatson: @quickactivate, datadir
@quickactivate "MCScatteringDataAnalysis"

using Dates: Dates, now
using Format: format

const seeds = 1:1000

const param_filename = "mc_in.txt"

const params = """
SKSPD    0 20  0                Shock speed as: km/sec, Lorentz factor, or beta.  1st nonzero value used
NITRS      20                   Number of iterations to perform
XNPER   1.0E2 1.0E4             Number of time steps per gyro period, for "coarse" and "fine" scattering
NIONS       2                   Number of different ion species.  1st MUST be protons.  Search "NIONS" for details.
    1     1  1E6  1               A(-99 for e-), Z, T0[K], & den0[/cc] for species 1
    -99   1  1E6  0               A(-99 for e-), Z, T0[K], & den0[/cc] for species 2
TELEC       0                   Far UpS elec temp.  Is set to 0 if e- is a separate species.
INDST       1                   Input distribution.  1 = thermal, 2 = delta function, 3 = other
ENINJ     1.0E3                 Injection energy[keV] and min PSD limit if "INDST" = 2.  Ignored otherwise.
INJWT       1                   How to assign initial weights.  1 = equal by particle, 2 = equal by bin
ENMAX    0 0 1.0E+10            Max particle energy as Emax[keV], Emax/nuc[keV/aa], pmax/(m_pc).  1st nonzero value used.  Also used in setting PSD
GYFAC       1                   Gyrofactor, i.e. MFP = gyrofac*gyroradius.  Ignored if "NWFRG" = 66
BMAGZ    1.0E-5                 Far upstream magnetic field in Gauss
THTBZ       0                   Upstream angle[deg] between B-field and shock normal
XGDUP      -1.0E7               Start of grid[rg0].  Can not be closer than UpS FEB.
XGDDW      10                   Downstream limit of grid[rg0]
FEBUP   -1.0E2  0               Upstream FEB, in [rg0] or [pc].  1st non-zero value used.  Must be negative.
FEBDW    0.0    0               Downstream FEB, in [rg0] or [pc].  1st non-zero value used.  Must be positive if included; set <= 0 otherwise.
NSPEC       0                   Number of x-positions [rg0] where particle spectrum is calculated
DNDPS      66                   Enter "66" to write a separate dNdp for each iteration
NWFRG       0                   Enter "66" to define custom f(r_g) in subr "scattering"; else = eta_mfp*r_g
NPTLO      600   1200           Number of particles to inject, and target number of particles at low-E pcuts.
NPTHI     2000 1.0E+06          Target number of particles at hi-E pcuts, and cutoff kinetic energy [keV/aa] between the two
PCUTS                           List of momentum cutoffs [aa*m_pc] to use during iterations.
    1.000E-02                     pcut  1
    6.000E-01                     pcut  2
    1.600E+00                     pcut  3
    2.000E+00                     pcut  4
    4.500E+00                     pcut  5
    9.000E+00                     pcut  6
    3.000E+01                     pcut  7
    5.000E+01                     pcut  8
    2.000E+02                     pcut  9
    3.000E+02                     pcut  10
    5.000E+02                     pcut  11
    1.000E+03                     pcut  12
    2.000E+03                     pcut  13
    5.000E+03                     pcut  14
    1.000E+04                     pcut  15
    3.162E+04                     pcut  16
    1.000E+05                     pcut  17
    3.162E+05                     pcut  18
    1.000E+06                     pcut  19
    3.162E+06                     pcut  20
    1.000E+07                     pcut  21
    1.778E+07                     pcut  22
    3.162E+07                     pcut  23
    5.623E+07                     pcut  24
    1.000E+08                     pcut  25
    1.778E+08                     pcut  26
    3.162E+08                     pcut  27
    5.623E+08                     pcut  28
    1.000E+09                     pcut  29
    1.778E+09                     pcut  30
    3.162E+09                     pcut  31
    5.623E+09                     pcut  32
    1.000E+10                     pcut  33
    1.778E+10                     pcut  34
    3.162E+10                     pcut  35
    5.623E+10                     pcut  36
    1.000E+11                     pcut  37
    1.778E+11                     pcut  38
    3.162E+11                     pcut  39
    5.623E+11                     pcut  40
    1.000E+12                     pcut  41
    1.778E+12                     pcut  42
    3.162E+12                     pcut  43
    5.623E+12                     pcut  44
    1.000E+13                     pcut  45
    -1                            pcut  46
NOSHK       0                   Enter "66" for no shock, i.e. force r_comp = 1.  For testing only.
NOSCT       0                   Enter "66" for scatter-FREE propagation.  For testing only.
NODSA       0                   Enter "66" to turn off DSA (diffusive shock acceleration).  For testing only.
SMSHK      66                   Enter "66" to keep velocity profile const. between iterations (i.e. NO smoothing)
SMIWT       1                   Factor to weight old profile with in shock smoothing. > 1 favors old profile over new.  = 1 to use mean of the two
SMVWT       0                   Enter "66" to increase old profile weighting for later iterations
SMMOE       0                   For profile iteration, = 0 uses only mom. flux eq, = 1 only energy.  Must be in [0,1]
SMPFP       0                   For profile iteration, = 0 uses only flux to find pressure, = 1 only PSD data.  Must be in [0,1]
RCOMP      -1                   Target compression ratio.  = -1 to use R-H value of 3.0583
OLDIN       0                   Enter "66" to read in old profile from "mc_grid_old.dat"
OLDDT    1300  5 100            If reading old profile, number of lines to skip, number of profiles to average, number of lines per profile
AGEMX   3.15E11                 Maximum allowed cosmic ray age (sec, explosion frame).  Ignored if <= 0
TCUTS                           List of time cutoffs to use for tracking particles
    1E03                          tcut  1
    1E04                          tcut  2
    1E05                          tcut  3
    1E06                          tcut  4
    1E07                          tcut  5
    1E08                          tcut  6
    1E09                          tcut  7
    1E10                          tcut  8
    1E11                          tcut  9
    3E13                          tcut  10
    -1                            tcut  11
RETRO      66                   Enter "66" to use retro time calc downstream.  Ignored if "AGEMX" < 0
FPUSH      66                   Enter "66" for fast upstream transport
FPSTP      -1                   UpS x position[rg0] where PROTON fast transport stops
ARTSM       0 0                 Artificial smoothing: start position[rg0], scale factor.  Ignored if 1st input >= 0
EMNFP   1.0E4                   Kinetic energy[keV] below which electrons have constant mfp.  Ignored if <= 0
NORAD       0                   Enter "66" for NO radiation losses
PHOTN       0                   Enter "66" to calculate photon production
JETRD   0.438                   Jet shock radius[pc] used for normalizing photon production.  Ignored if "PHOTN" != 66
JETFR       0 5                 Frac of sphere producing CRs, or jet open ang[deg].  1st non-zero value used.  Ignored if "PHOTN" != 66
JETDS   1.0E6                   Distance from jet to observer[kpc].  Used to calc photon fluxes & CMB parameters
ENXFR     0.1                   Fraction of ion energy transferred to electrons at 1st shock crossing. Must be [0,1]; 0 to ignore
NSHLS       5 2                 Number of upstream & downstream shells for photon emission.  Ignored if "PHOTN" != 66
BTRBF       1                   Fraction of compressed B turb to use in scattering & losses.  Must be in [0,1].  0 means ignore
BAMPF       1                   Amplification factor to use for B field.  1.0 means no amp
NWEPB      66                   Enter "66" to define custom epsilon_B on grid; else = B(x) controlled by "BTRBF" and "BAMPF"
PSDBD    10  10                 Number of PSD bins per decade in (1) momentum, (2) theta
PSDTB   119   4                 Number of PSD ang bins: (1) # linear (cosine) bins, (2) # log (theta) decades
ENDIN
 **** End of input data ****
"""

const runpath = datadir("Lorentz-5-raw")

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
    return
end

function initdirs(seeds)
    mkpath(runpath)
    cd(runpath)
    for seed in seeds
        dirname = format("Seed-{:0>4}", seed)
        # Create the folder
        mkdir(dirname)

        # `cd()` into the folder and write the parameters
        cd(dirname) do
            writeparams(param_filename, seed)
            @info("Created parameter file for seed = $seed at $(pwd())")
        end
    end
    return
end

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
                stdout = "logs-stdout.txt", stderr = "logs-stderr.txt"
            ) |> run
        end
        time = Dates.format(now(), timefmt)
        @info("Ended run in $dirname ($time)")
    end
    return
end

end
