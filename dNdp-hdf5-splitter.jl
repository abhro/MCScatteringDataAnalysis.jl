using JLD2
using DataFrames

"""
dehistogram each dataframe of a GroupedDataFrame, then return one big DataFrame
containing all the rows.
"""
function dehistogram(gdf::GroupedDataFrame)
    df_dehistogrammed = DataFrame()
    for df in gdf
        # de-histogram each of the sub-dataframes (make a copy first)
        df_dehist = df[begin+1:2:end,:]

        # XXX XXX XXX XXX
        # **Need to explore `psd_mom_bounds`**, and how data repeats in adjacent cells (high/low bounds)
        # -> Choose the upper bounds for de-histogramming the data, because the lowest
        #    bin covers -99.0 to -19.3 g⋅cm/s, which is too wide. The highest bin covers
        #    -2.3997 to -2.2997 g⋅cm/s
        # -> runs into a fencpost problem. complete last cell (343 for each grouped df)
        #    has a unrepeated value?

        # push that copy into a bigger temp df
        # the cols attribute is an additional sanity check
        append!(df_dehistogrammed, df_dehist, cols = :orderequal)
    end
    return df_dehistogrammed
end

CR_df = let
    df = load_object("dNdp-CR.jld2")
    gdf = groupby(df, [:initial_seed, :iteration])
    insertcols!(df, 1, :run_id => groupindices(gdf))
    # drop the two columns which are now redundant (or not ig)
    #select!(df, Not(:initial_seed, :iteration))
    # column i is pointless, only has 66s
    # plot is useful for pgf plot, not here
    select!(df, Not(:i, :plot))

    disallowmissing!(df, error = false)

    df
end

CR_p_df, CR_e_df = let
    # split CR_df based on ion and iteration number
    # split therm_df based on ion and iteration number because each of them is
    # a complete run which needs to be de-histogrammed
    gdf = groupby(CR_df, [:iter, :ion])
    CR_df_dehistogrammed = dehistogram(gdf)

    # if all three dNdp_cr is missing, drop that row
    subset!(
        CR_df_dehistogrammed,
        [:log_dNdp_cr_sf, :log_dNdp_cr_pf, :log_dNdp_cr_ISM] =>
            ByRow((cols...) -> !all(ismissing, cols)),
    )

    # group that bigger temp df by ions and return it
    groupby(CR_df_dehistogrammed, :ion)
end

# Separate each of the proton and electron DataFrames by momentum.
CR_p_gdf_momentum = groupby(CR_p_df, :log_p_nat)
CR_e_gdf_momentum = groupby(CR_e_df, :log_p_nat)
save_object("dNdp-CR-protons-momentum-split.jld2", CR_p_gdf_momentum)
@info "Saved dNdp-CR-protons-momentum-split.jld2"
save_object("dNdp-CR-electrons-momentum-split.jld2", CR_e_gdf_momentum)
@info "Saved dNdp-CR-electrons-momentum-split.jld2"

# Separate each of the proton and electron DataFrames based on iteration.
CR_p_gdf_iteration = groupby(CR_p_df, :iter)
CR_e_gdf_iteration = groupby(CR_e_df, :iter)
save_object("dNdp-CR-protons-iteration-split.jld2", CR_e_gdf_iteration)
@info "Saved dNdp-CR-protons-iteration-split.jld2"
save_object("dNdp-CR-electrons-iteration-split.jld2", CR_e_gdf_iteration)
@info "Saved dNdp-CR-electrons-iteration-split.jld2"


# repeat for thermal dNdp
#therm_df = let
#    df = load_object("dNdp_therm.jld2")
#    gdf = groupby(df, [:initial_seed, :iteration])
#    insertcols!(df, 1, :iter => groupindices(gdf))
#    # drop the two columns which are now redundant (or not ig)
#    #select!(df, Not(:initial_seed, :iteration))
#    # column i is pointless, only has 66s
#    # plot is useful for pgf plot, not here
#    select!(df, Not(:i, :plot))
#    disallowmissing!(df, error = false)
#    df
#end

#therm_p_df, therm_e_df = let
#    # split therm_df based on ion and iteration number because each of them is
#    # a complete run which needs to be de-histogrammed
#    gdf = groupby(therm_df, [:iter, :ion])
#    therm_df_dehistogrammed = dehistogram(gdf)
#
#    # if all three dNdp_therm is missing, drop that row
#    subset!(
#        therm_df_dehistogrammed,
#        [:log_dNdp_therm_sf, :log_dNdp_therm_pf, :log_dNdp_therm_ISM] =>
#            ByRow((cols...) -> !all(ismissing, cols)),
#    )

#    # group that bigger temp df by ions and return it
#    groupby(therm_df_dehistogrammed, :ion)
#end

## Separate each of the proton and electron DataFrames by momentum.
## XXX there's not a single specific column of momentum to split it on
##therm_p_gdf_momentum = groupby(therm_p_df, :log_p_nat)
##therm_e_gdf_momentum = groupby(therm_e_df, :log_p_nat)
##save_object("dNdp-therm-protons-momentum-split.jld2", therm_p_gdf_momentum)
##save_object("dNdp-therm-electrons-momentum-split.jld2", therm_e_gdf_momentum)

## Separate each of the proton and electron DataFrames based on iteration.
#therm_p_gdf_iteration = groupby(therm_p_df, :iter)
#therm_e_gdf_iteration = groupby(therm_e_df, :iter)
#save_object("dNdp-therm-protons-iteration-split.jld2", therm_p_gdf_iteration)
#save_object("dNdp-therm-electrons-iteration-split.jld2", therm_e_gdf_iteration)
