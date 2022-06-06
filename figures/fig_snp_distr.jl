using DataFrames,CSV,StatsBase,RCall

snp = CSV.read("/Users/aaroneisman/master_snp_maf_union_gwas_df_v4.csv",DataFrame)
snp_jhs = filter(row -> row.jhs_select == 1, snp)
snp_mesa = filter(row -> row.mesa_select == 1, snp)
snp_her = filter(row -> row.her_select == 1, snp)

jhs_counts = collect(values(countmap(snp_jhs.prot)))
mesa_counts = collect(values(countmap(snp_mesa.prot)))
her_counts = collect(values(countmap(snp_her.prot)))

@rput(jhs_counts,mesa_counts,her_counts)
R"""
    fn = "jhs_hist_counts.pdf"
    pdf(fn, width = 7, height = 7)
    hist(jhs_counts, main = "JHS IVs / Protein", breaks = 0:10)
    dev.off()


    fn = "mesa_hist_counts.pdf"
    pdf(fn, width = 7, height = 7)
    hist(mesa_counts, main = "MESA IVs / Protein", breaks = 0:10)
    dev.off()

    fn = "her_hist_counts.pdf"
    pdf(fn, width = 7, height = 7)
    hist(her_counts, main = "HERITAGE IVs / Protein", breaks = 0:10)
    dev.off()
"""