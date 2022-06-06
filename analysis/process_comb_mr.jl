using DataFrames,CSV,Dates

header_names = 
["mr_file","mr_folder","mr_str","B_IVW","SE_IVW","CI95_L_IVW","CI95_U_IVW","P_IVW","B_MaxLik","SE_MaxLik","CI95_L_MaxLik","CI95_U_MaxLik","P_MaxLik","B_Simple_median","SE_Simple_median","CI95_L_Simple_median","CI95_U_Simple_median","P_Simple_median","B_Weighted_median","SE_Weighted_median","CI95_L_Weighted_median","CI95_U_Weighted_median","P_Weighted_median","B_MRE","SE_MRE","CI95_L_MRE","CI95_U_MRE","P_MRE","B_MRE_intercept","SE_MRE_intercept","CI95_L_MRE_intercept","CI95_U_MRE_intercept","P_MRE_intercept"]

files_arr = ["mr.comb.gwas_protein_out_heritage_fdr_bonf_ld_0.001_tss_1000000_maf_0.01.csv","mr.comb.gwas_protein_out_jhs_fdr_bonf_ld_0.001_tss_1000000_maf_0.01.csv","mr.comb.gwas_protein_out_mesa_fdr_bonf_ld_0.001_tss_1000000_maf_0.01.csv"]

for file in files_arr

    mr = CSV.read(file,DataFrame, header = header_names)

    filter!(row -> !ismissing(row.mr_str),mr)

    mr[!,:exp] .= ""
    mr[!,:out] .= ""

    for i in 1:size(mr)[1]
        t_mr_str_arr = split(mr.mr_str[i],"->")
        t_exp = t_mr_str_arr[1]
        t_out = t_mr_str_arr[2]

        mr.exp[i] = t_exp
        mr.out[i] = t_out
    end

    CSV.write(string(today())*"."*file,mr)
end

#    filter(row -> row.mr_str == "SL005574->MN1074", mr)
#    filter(row -> row.mr_str == "SL000668->MP1107", mr)
