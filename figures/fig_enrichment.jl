using DataFrames,CSV,RCall
try
    include("../analysis/omics_tools.jl")
    include("../analysis/enrichment_tools.jl")
catch
    include("./analysis/omics_tools.jl")
    include("./analysis/enrichment_tools.jl")
end

#### Enrichment Figure Panels

conf = "./figures/fig.conf";
project_path = read_conf(conf,"project_path");
results_sub_path = read_conf(conf,"results_sub_path");
annot_sub_path = read_conf(conf,"annot_sub_path");
fig_enrich_out = read_conf(conf,"fig_enrich_out");

ogm = CSV.read(project_path*annot_sub_path*"ogm_meta_key_anno.csv",DataFrame);
pm_class_set = Set(ogm.pm_class);
pm_subclass_set = Set(ogm.pm_subclass);

for class in pm_class_set
    ogm[!,Symbol(class)] = ogm.pm_class .== class
end

for class in pm_subclass_set
    ogm[!,Symbol(class)] = ogm.pm_subclass .== class
end

## A) APOE correlation histogram
#cor_qval_df = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/data/multiomics/v4b_jhs_mesa_her_cor_filt_qval_by_protein.csv",DataFrame)
#cor_qval_df = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/data/multiomics/jhs_mesa_her_cor_filt_qval_by_protein.csv",DataFrame)
cor_df = CSV.read(project_path*results_sub_path*"all_cor_filt_meta_qval.csv",DataFrame);
filter!(row -> !ismissing(row.meta_p_qval),cor_df);
for class in union(pm_class_set,pm_subclass_set)
#for class in pm_class_set
    cor_df[!,Symbol(class)] .= 0
    t_class_set = Set(ogm.ogm_id[Bool.( ogm[:,Symbol(class)] )])
    for i in 1:size(cor_df)[1]
        if in(cor_df.OGM[i],t_class_set)
            cor_df[i,Symbol(class)] = 1
        end
    end
end

apoe_cor_df = filter(row -> row.PROT == "SL000276",cor_df)
apoe_all_rho = apoe_cor_df.meta_rho
sig = 0.05
#apoe_lipid_cor = filter(row -> row.meta_p_qval < sig, apoe_cor_df).meta_rho
apoe_lipid_rho = filter(row -> row.Combined_Lipids == 1, apoe_cor_df).meta_rho
apoe_lipid_rho_sig = filter(row -> row.Combined_Lipids == 1 && row.meta_p_qval < sig, apoe_cor_df).meta_rho

function volcano_plot(df,p_col,b_col,sig_var,fn)
    #sort!(df,Symbol(sig_var))
    df = filter(row -> !ismissing(row[Symbol(p_col)]), df)
    df[!,:log10_p] = -log.(10,df[:,Symbol(p_col)])

    sig = 0.05
    df_not_sig = filter(row -> row[Symbol(p_col)] > sig, df)
    df_lipid_sig = filter(row -> row[Symbol(p_col)] < sig && row[Symbol(sig_var)] == 1, df)
    df_other_sig = filter(row -> row[Symbol(p_col)] < sig && row[Symbol(sig_var)] == 0, df)

    log10_p = -log.(10,df[:,Symbol(p_col)])
    betas = df[:,Symbol(b_col)]
    sig = (df[:,Symbol(sig_var)] .== 1) .+ 1
    v_plot_data = hcat(betas,log10_p)

    v_plot_data = DataFrame(betas = betas, log10_p = log10_p)

    #v_plot_data_not_sig = filter(row -> row.color == "gray50")
    @rput(v_plot_data,df_not_sig, df_lipid_sig, df_other_sig, fn)
    R"""
        pdf(fn, width = 4.75, height = 4.75)

        #plot(v_plot_data$log10_p~v_plot_data$betas,xlim = c(-max(abs(v_plot_data$betas)),max(abs(v_plot_data$betas))), pch = 6, xlab = "META-Analyzed Beta Coefficient", ylab = "-log10(q-value)")
        plot(df_not_sig$log10_p~df_not_sig$meta_rho,xlim = c(-.5,.5), ylim = c(0,300), pch = 1, xlab = "APOE:Metabolite Correlation Coefficient", ylab = "-log10(q-value)", col = "gray50")
        
        points(df_other_sig$meta_rho,df_other_sig$log10_p)
        points(df_lipid_sig$meta_rho,df_lipid_sig$log10_p, pch = 19, col = "green3")
        abline(h = 1.3, lty = 2, lwd = 2, col = "gray50")
        legend(x = "topleft", legend = (c("APOE","Lipid Associations","Other Associations", "Not Significant (q > 0.05)")), col = c("white","green3","black","gray50"), pch = c(1,19,1,1), bty = "n")
        #points(v_plot_data$log10_p~v_plot_data$betas, col = v_plot_data$color, pch = v_plot_data$shape)
        

        dev.off() 
    """
end
df = apoe_cor_df
p_col = "meta_p_qval"
b_col = "meta_rho"
sig_var = "Combined_Lipids"
fn = project_path*fig_enrich_out*"apoe_lipid_volcano.pdf"
volcano_plot(apoe_cor_df, p_col, b_col, sig_var, fn)

#=@rput(apoe_all_rho,apoe_lipid_rho,apoe_lipid_rho_sig)
R"""
    b = c(-1,-6:6*.05,1)
    b = c(-1,-3:3*.1,1)
    b = 30

    fn = "apoe_rho_hist.pdf"
    pdf(fn, width = 6, height = 6)
    
    hist(apoe_all_rho, col = "white", xlim = c(-0.6,0.6), breaks = b, xlab = "ApoE:Metabolite Correlation Coefficient", main = "")
    hist(apoe_lipid_rho, col = "green3", add = TRUE, breaks = b, density = 30)
    hist(apoe_lipid_rho_sig, col = "green3", add = TRUE, breaks = b)

    legend(x = "topright",legend = c("Lipids (q < 0.05)","Lipids (q > 0.05)","Other Metabolites"), fill = c("green3","green3","white"), density = c(100,30,100))

    dev.off()

"""
=#

## B) APOE Enrichment Methods
 apoe_cor_lipid_hits = sort(apoe_cor_df,:meta_p).Combined_Lipids
 fn = project_path*fig_enrich_out*"apoe_lipid_hits.pdf"
 @rput(apoe_cor_lipid_hits,fn)
 R"""
    pdf(fn, width = 6, height = 6)
    plot(1, type="n", xlab="", ylab="", xlim = c(0,400), ylim=c(0, 10))
    for (i in 1:length(apoe_cor_lipid_hits)){
        if (apoe_cor_lipid_hits[i] == 1){
            abline(v = i, col = "green3")
            #abline(v = i, col = "#377c78")
        }
        abline(v = length(apoe_cor_lipid_hits)+5, col = "red")
    }
    dev.off()
 """

sort!(apoe_cor_df,:meta_p)
apoe_arr_for_ogea = Matrix(apoe_cor_df[:,[Symbol(:Combined_Lipids),:meta_rho]])
p = 1
apoe_var_cum_sum = ogea(apoe_arr_for_ogea,p)
apoe_ogea_graph_df = DataFrame(ApoE = apoe_var_cum_sum)
fn = project_path*fig_enrich_out*"apoe_ogea_graph.pdf"
ogea_graph(apoe_ogea_graph_df, ["ApoE"], fn, ["ApoE"], ["black"],0)

#=ogm.var2_N_Acetyl_AA = contains.(ogm.ogm_name,r"N-[Aa]cetyl")
var2_N_Acetyl_AA_set = Set(filter(row -> row.var2_N_Acetyl_AA, ogm).ogm_id)

acy1_cor_df = filter(row -> row.PROT == "SL005574",cor_qval_df)
acy1_cor_df[!,:var2_N_Acetyl_AA] .= 0
for i in 1:size(acy1_cor_df)[1]
    if in(acy1_cor_df.jhs_var2_alt[i],var2_N_Acetyl_AA_set)
        acy1_cor_df.var2_N_Acetyl_AA[i] = 1
    end
end

acy1_arr_for_ogea = Matrix(acy1_cor_df[:,[Symbol(:var2_amino_acid),:meta_rho]])
acy1_var_cum_sum = ogea(acy1_arr_for_ogea,p)
acy1_ogea_graph_df = DataFrame(Acy1 = acy1_var_cum_sum)
ogea_graph(acy1_ogea_graph_df, ["Acy1"], "acy1_ogea_graph.pdf", ["Acy1"], ["black"],0)
=#