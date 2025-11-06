using StaticArrays: SVector

#=
This file lays out the columns and data types for various fixed width format
files written by the mc_cr program
=#

# from Fortran subroutine `print_dNdp_esc`
# the Fortran file unit in use is 527
const esc_cols = SVector(
    ColumnSpecification(:itr, Int),
    ColumnSpecification(:plot, Int),
    ColumnSpecification(:ion, Int),
    ColumnSpecification(:log_p_cgs, Float64),
    ColumnSpecification(:log_p_nat, Float64),
    ColumnSpecification(:log_dNdp_esc_UpS_SF, Float64, true),
    ColumnSpecification(:log_dNdp_esc_UpS_PF, Float64, true),
    ColumnSpecification(:log_dNdp_esc_UpS_IF, Float64, true),
    ColumnSpecification(:log_dNdp_esc_DwS_SF, Float64, true),
    ColumnSpecification(:log_dNdp_esc_DwS_PF, Float64, true),
    ColumnSpecification(:log_dNdp_esc_DwS_IF, Float64, true),
)

# from Fortran subroutine `smooth_grid_par`, 312 lines after the function signature
# the Fortran file unit in use is 101
const grid_cols = SVector(
    ColumnSpecification(:itr, Int),
    ColumnSpecification(:grid_index, Int), # index of grid zone, called `i` in Fortran code

    ColumnSpecification(:x_grid_rg, Float64),
    ColumnSpecification(:x_grid_log, Float64),
    ColumnSpecification(:x_grid_cm, Float64),
    ColumnSpecification(:x_grid_log_cm, Float64),

    # momentum flux
    ColumnSpecification(:pxx_norm, Float64),
    ColumnSpecification(:pxx_norm_log, Float64),
    ColumnSpecification(:pxz_norm, Float64),
    ColumnSpecification(:pxz_norm_log, Float64),

    # energy flux
    ColumnSpecification(:en_norm, Float64),
    ColumnSpecification(:en_norm_log, Float64),

    # velocity normalized to first grid zone
    ColumnSpecification(:ux_norm, Float64),
    ColumnSpecification(:ux_norm_log, Float64),
    ColumnSpecification(:uz_norm, Float64),
    ColumnSpecification(:uz_norm_log, Float64),

    ColumnSpecification(:bmag, Float64),
    ColumnSpecification(:bmag_log, Float64),

    ColumnSpecification(:θ_deg, Float64), # angle between magnetic field and shock normal

    ColumnSpecification(:γᵤ_sf, Float64), # lorentz factor at shock frame

    # density ratio = (γ_Z β_Z) / (γ β)
    ColumnSpecification(:inv_density_ratio, Float64),
    ColumnSpecification(:density_ratio, Float64),

    ColumnSpecification(:pressure_px_log, Float64),
    ColumnSpecification(:pressure_en_log, Float64),
    ColumnSpecification(:pressure_psd_par_log, Float64),
    ColumnSpecification(:pressure_psd_perp_log, Float64),
    ColumnSpecification(:pressure_tot_MC_log, Float64),
    ColumnSpecification(:pressure_aniso, Float64),
    # tp -> test particle limit, absence of DSA
    ColumnSpecification(:pressure_px_tp, Float64),
    ColumnSpecification(:pressure_en_tp, Float64),
    ColumnSpecification(:pressure_Z, Float64),

    ColumnSpecification(:log_q_remaining_cal_px, Float64),
    ColumnSpecification(:log_q_remaining_cal_en, Float64),
    ColumnSpecification(:εB_grid, Float64),
    ColumnSpecification(:log_εB_grid, Float64),
)

# from Fortran subroutine `get_normalized_dNdp`, 274 lines after the function signature
# the Fortran file unit in use is 507
const therm_cols = SVector(
    ColumnSpecification(:i, Int),
    ColumnSpecification(:plot, Int),
    ColumnSpecification(:ion, Int),
    ColumnSpecification(:log_dNdp_therm_pvals_cgs_sf, Float64, true),
    ColumnSpecification(:log_dNdp_therm_pvals_nat_sf, Float64, true),
    ColumnSpecification(:log_dNdp_therm_sf, Float64, true),
    ColumnSpecification(:log_dNdp_therm_pvals_cgs_pf, Float64, true),
    ColumnSpecification(:log_dNdp_therm_pvals_nat_pf, Float64, true),
    ColumnSpecification(:log_dNdp_therm_pf, Float64, true),
    ColumnSpecification(:log_dNdp_therm_pvals_cgs_ISM, Float64, true),
    ColumnSpecification(:log_dNdp_therm_pvals_nat_ISM, Float64, true),
    ColumnSpecification(:log_dNdp_therm_ISM, Float64, true),
)

# from Fortran subroutine `get_normalized_dNdp`, 343 lines after the function signature
# the Fortran file unit in use is 517
const CR_cols = SVector(
    ColumnSpecification(:i, Int),
    ColumnSpecification(:plot, Int),
    ColumnSpecification(:ion, Int),
    ColumnSpecification(:log_p_cgs, Float64),
    ColumnSpecification(:log_p_nat, Float64),
    ColumnSpecification(:log_dNdp_cr_sf, Float64, true),
    ColumnSpecification(:log_dNdp_cr_pf, Float64, true),
    ColumnSpecification(:log_dNdp_cr_ISM, Float64, true),
    ColumnSpecification(:log_dNdp_all_sf, Float64, true),
    ColumnSpecification(:log_dNdp_all_pf, Float64, true),
    ColumnSpecification(:log_dNdp_all_ISM, Float64, true),
)


# from Fortran subroutine `tcut_print`
# the Fortran file unit in use is 670
const coupled_weights_cols = SVector(
    ColumnSpecification(:itr, Int),
    ColumnSpecification(:plot, Int),
    ColumnSpecification(:log_tcut, Float64),
    # the following number of columns are based on the number of species in mc_in.txt.
    # right now, we are concerned with two: protons and electrons. the column
    # specification will need to be changed if mc_in.txt is changed
    ColumnSpecification(:log_weight_proton, Float64),
    ColumnSpecification(:log_weight_electron, Float64),
)

const n_tcuts = 10
# from Fortran subroutine `tcut_print`
# the Fortran file unit in use is 671
const coupled_spectra_cols = SVector(
    ColumnSpecification(:itr, Int),
    ColumnSpecification(:plot, Int),
    ColumnSpecification(:ion, Int),

    ColumnSpecification(:log_p_cgs, Float64),
    ColumnSpecification(:log_p_nat, Float64),

    # the following number of columns are based on the number of tcuts in mc_in.txt.
    ColumnSpecification.(
        Symbol.("log_spectra_tcut_", Base.OneTo(n_tcuts)),
        Float64,
        true,
    )...,
)
