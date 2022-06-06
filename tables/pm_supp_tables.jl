using CSV, DataFrames, RCall, Distributions, StatsBase, Printf
try
    include("../analysis/omics_tools.jl")
catch
    include("./analysis/omics_tools.jl")
end

conf = "./tables/tab.conf";
project_path = read_conf(conf,"project_path");
data_sub_path = read_conf(conf,"data_sub_path");
tab_out_sub_path = read_conf(conf,"tab_out_sub_path");
results_sub_path = read_conf(conf,"results_sub_path");


### Supplementary Table 1 - patient characteristics
# n / age / sex/ bmi / egfr / ethnicity
function mean_std_str(df,var)
    t_mean = @sprintf("%.1f",mean(skipmissing(df[:,Symbol(var)])))
    t_sd = @sprintf("%.1f",std(skipmissing(df[:,Symbol(var)])))
    mean_std_str = "$(t_mean) ($(t_sd))"

    return mean_std_str
end

function n_pct_str(df,var,ref_val)
    t_dict = countmap(df[:,Symbol(var)])
    t_n = t_dict[ref_val]
    t_pct = t_n / sum(values(t_dict)) * 100

    t_pct_str = @sprintf("%.0f",t_pct)

    n_pct_str = "$(t_n) ($(t_pct_str)%)" 

    return n_pct_str
end

s1 = DataFrame(Characteristic = ["N","Age, years","Female","BMI, kg/m2","eGFR, ml/min/1.73m2","Race/Ethnicity","Black","Asian","Hispanic","White"])

# JHS
jhs_adj_file = read_conf(conf,"jhs_adj_file");
jhs = CSV.read(project_path*data_sub_path*jhs_adj_file,DataFrame);
jhs_incl_bin = completecases(jhs[:,[:age,:sex,:SL005574,Symbol("C18:1 SM")]])
jhs_incl = jhs[jhs_incl_bin,[:age,:sex,:bmi,:egfr]]

jhs_n_str = string(size(jhs_incl)[1])
jhs_age_str = mean_std_str(jhs_incl,"age")
jhs_n_female = sum(jhs_incl.sex .== "Female")
jhs_female_str = n_pct_str(jhs_incl,"sex","Female")
jhs_bmi_str = mean_std_str(jhs_incl,"bmi")
jhs_egfr_str = mean_std_str(jhs_incl,"egfr")
jhs_black_str = "$(jhs_n_str) (100%)"

s1[!,:JHS] = [jhs_n_str,jhs_age_str,jhs_female_str,jhs_bmi_str,jhs_egfr_str,"",jhs_black_str,"-","-","-"]

# MESA
mesa_adj_file = read_conf(conf,"mesa_adj_file")
mesa = CSV.read(project_path*data_sub_path*mesa_adj_file,DataFrame);
mesa_incl = filter(row -> !ismissing(row.PlateId), mesa)
mesa_n_str = string(size(mesa_incl)[1])
mesa_age_str = mean_std_str(mesa_incl,"age")
mesa_female_str = n_pct_str(mesa_incl,"sex",0) #0=FEMALE. 1=MALE
mesa_bmi_str = mean_std_str(mesa_incl,"bmi")
mesa_egfr_str = mean_std_str(mesa_incl,"egfr")

# 1=WHITE, CAUCASIAN 2=CHINESE AMERICAN 3=BLACK, AFRICAN-AMERICAN 4=HISPANIC
mesa_black_str = n_pct_str(mesa_incl,"race1c",3)
mesa_chinese_str = n_pct_str(mesa_incl,"race1c",2)
mesa_hispanic_str = n_pct_str(mesa_incl,"race1c",4)
mesa_white_str = n_pct_str(mesa_incl,"race1c",1)

s1[!,:MESA] = [mesa_n_str,mesa_age_str,mesa_female_str,mesa_bmi_str,mesa_egfr_str,"",mesa_black_str,mesa_chinese_str,mesa_hispanic_str,mesa_white_str]

# HERITAGE

her_adj_file = read_conf(conf,"her_adj_file")
her = CSV.read(project_path*data_sub_path*her_adj_file,DataFrame);
her_incl = her
her_n_str = string(size(her_incl)[1])
her_age_str = mean_std_str(her_incl,"age")
her_female_str = n_pct_str(her_incl,"sex",2) #1=M, 2=F
her_bmi_str = mean_std_str(her_incl,"bmi")
her_egfr_str = "-"

# RACE, 1=Black;2=Caucasian
her_black_str = n_pct_str(her_incl,"RACE",1)
her_chinese_str = "-"
her_hispanic_str = "-"
her_white_str = n_pct_str(her_incl,"RACE",2)

s1[!,:HERITAGE] = [her_n_str,her_age_str,her_female_str,her_bmi_str,her_egfr_str,"",her_black_str,her_chinese_str,her_hispanic_str,her_white_str]

CSV.write(project_path*tab_out_sub_path*"s1_characteristics.csv",s1)


### Supplementary Table 8
### CHANGE TO entrez_gene_symbol, REMOVE -### in JHS and MESA
jhs_prot_file = read_conf(conf,"jhs_prot_file");
mesa_prot_file = read_conf(conf,"mesa_prot_file");
her_prot_file = read_conf(conf,"her_prot_file");
her_conv_file = read_conf(conf,"her_conv_file");
cor_data_file = read_conf(conf,"cor_data_file");
kegg_annotations_file = read_conf(conf,"kegg_annotations_file");

cor_df = CSV.read(project_path*results_sub_path*cor_data_file,DataFrame)
prots_set = Set(filter(row -> !ismissing(row.meta_p), cor_df).PROT)

kegg_annot = CSV.read(project_path*data_sub_path*kegg_annotations_file,DataFrame)
kegg_annot = rename(kegg_annot,Dict(:entrez_gene_symbol => "Target_Gene", :target_full_name => "Target_Name", :soma_id => "SOMA_ID"))

jhs_prot = CSV.read(project_path*data_sub_path*jhs_prot_file,DataFrame)
mesa_prot = CSV.read(project_path*data_sub_path*mesa_prot_file,DataFrame)
her_prot = CSV.read(project_path*data_sub_path*her_prot_file,DataFrame)
her_conv = CSV.read(project_path*data_sub_path*her_conv_file,DataFrame)

sort!(filter!(row -> in(row.SomaId,prots_set), jhs_prot),:SomaId)
select!(jhs_prot,[:SomaId,:SeqId,:TargetFullName])
jhs_prot = rename(jhs_prot,[:SOMA_ID,:JHS_Somalogic_Aptamer_ID,:JHS_Target_Name])

select!(mesa_prot,[:SomaId,:SeqId,:TargetFullName])
mesa_prot = rename(mesa_prot,[:SOMA_ID,:MESA_Somalogic_Aptamer_ID,:MESA_Target_Name])

her_sasname_keep_set = Set(her_conv.heritage_name)
filter!(row -> startswith(row.SomaId,"SL"),her_prot)
filter!(row -> in(row.SASNAME,prots_set) || in(row.SASNAME,her_sasname_keep_set),her_prot)
select!(her_prot,[:SomaId,:SeqId,:TargetFullName])
her_prot = rename(her_prot,[:SOMA_ID,:HERITAGE_Somalogic_Aptamer_ID,:HERITAGE_Target_Name])

jhs_mesa_prot = leftjoin(jhs_prot,mesa_prot, on = :SOMA_ID)

jhs_mesa_her_prot = leftjoin(jhs_mesa_prot,her_prot, on = :SOMA_ID)

jhs_mesa_her_prot = leftjoin(jhs_mesa_her_prot,select(kegg_annot,[:SOMA_ID,:Target_Gene,:Target_Name]), on = :SOMA_ID)

for i in 1:size(jhs_mesa_her_prot)[1]
    t_jhs_apt = split(jhs_mesa_her_prot.JHS_Somalogic_Aptamer_ID[i],"_")[1]
    jhs_mesa_her_prot.JHS_Somalogic_Aptamer_ID[i] = t_jhs_apt

    t_mesa_apt = split(jhs_mesa_her_prot.MESA_Somalogic_Aptamer_ID[i],"_")[1]
    jhs_mesa_her_prot.MESA_Somalogic_Aptamer_ID[i] = t_mesa_apt
end

s8 = sort(select!(jhs_mesa_her_prot,[:SOMA_ID,:Target_Gene,:Target_Name,:JHS_Somalogic_Aptamer_ID,:MESA_Somalogic_Aptamer_ID,:HERITAGE_Somalogic_Aptamer_ID]),:SOMA_ID)
s8 = rename(s8,Dict(:SOMA_ID => "Somalogic_ID"))

# remove duplicate HERITAGE somalogic aptamer 8409-3B
filter!(row -> ismissing(row.HERITAGE_Somalogic_Aptamer_ID) || row.HERITAGE_Somalogic_Aptamer_ID != "8409-3B", s8)

CSV.write(project_path*tab_out_sub_path*"s8_proteins.csv",s8)



### Supplementary Table 9
mets_set = Set(cor_df.OGM)

# update ogm table to include additional refmet names
# load refmet update table
ogm_refmet_update_file = read_conf(conf,"ogm_refmet_update_file")
ogm_meta_file = read_conf(conf,"ogm_meta_file")

refmet_update = CSV.read(project_path*data_sub_path*ogm_refmet_update_file,DataFrame)
ogm_meta = CSV.read(project_path*data_sub_path*ogm_meta_file,DataFrame)

refmet_update_ogm_ids = Set(refmet_update.ogm_id)

for i in 1:size(ogm_meta)[1]
    t_ogm_id = ogm_meta.ogm_id[i]
    if in(t_ogm_id,refmet_update_ogm_ids)
        t_refmet = filter(row -> row.ogm_id == t_ogm_id, refmet_update).refmet_name[1]
        ogm_meta.refmet_name[i] = t_refmet
        println("updated refmet name")
    end
end

filter!(row -> in(row.ogm_id,mets_set), ogm_meta)
ogm_meta[!,:Metabolite_Name] = ogm_meta.refmet_name
refmet_missing = ismissing.(ogm_meta.refmet_name)
ogm_meta[refmet_missing,:Metabolite_Name] = ogm_meta.ogm_name[refmet_missing]

ogm_meta.Metabolite_Name = lowercase.(ogm_meta.Metabolite_Name)

CSV.write(project_path*tab_out_sub_path*"s9_metabolites.csv",sort(select(ogm_meta,:Metabolite_Name)))

### Supplementary Table 2 - meta-analyzed correlation data
cor_data_file = read_conf(conf,"cor_data_file");
cor_df = CSV.read(project_path*results_sub_path*cor_data_file,DataFrame)

prot_id_name = select()
cor_df = rename(cor_df,Dict(:PROT => "SOMA_ID", :OGM => "ogm_id"))

cor_df = leftjoin(cor_df,select(kegg_annot,[:SOMA_ID,:Target_Gene,:Target_Name]), on = :SOMA_ID)
cor_df = leftjoin(cor_df,select(ogm_meta,[:ogm_id,:Metabolite_Name]), on = :ogm_id)

cor_df = rename(cor_df,Dict(:SOMA_ID => "Somalogic_ID"))


s2 = select(cor_df,[:Somalogic_ID,:Target_Gene,:Target_Name,:Metabolite_Name,:meta_rho,:meta_p,:meta_se,:meta_p_qval,:meta_rho_bmi,:meta_p_bmi,:meta_se_bmi,:meta_p_bmi_qval,:meta_rho_bmi_egfr,:meta_p_bmi_egfr,:meta_se_bmi_egfr,:meta_p_bmi_egfr_qval,:jhs_rho,:jhs_p,:jhs_n,:jhs_p_qval,:jhs_rho_bmi,:jhs_p_bmi,:jhs_n_bmi,:jhs_p_bmi_qval,:jhs_rho_bmi_egfr,:jhs_p_bmi_egfr,:jhs_n_bmi_egfr,:jhs_p_bmi_egfr_qval,:mesa_rho,:mesa_p,:mesa_n,:mesa_p_qval,:mesa_rho_bmi,:mesa_p_bmi,:mesa_n_bmi,:mesa_p_bmi_qval,:mesa_rho_bmi_egfr,:mesa_p_bmi_egfr,:mesa_n_bmi_egfr,:mesa_p_bmi_egfr_qval,:her_rho,:her_p,:her_n,:her_p_qval,:her_rho_bmi,:her_p_bmi,:her_n_bmi,:her_p_bmi_qval])

s2_names_dict = Dict(
:meta_rho => "META_Rho",
:meta_p => "META_P",
:meta_se => "META_SE",
:meta_p_qval => "META_Q",
:meta_rho_bmi => "META_Rho_BMI",
:meta_p_bmi => "META_P_BMI",
:meta_se_bmi => "META_SE_BMI",
:meta_p_bmi_qval => "META_Q_BMI",
:meta_rho_bmi_egfr => "META_Rho_BMI_eGFR",
:meta_p_bmi_egfr => "META_P_BMI_eGFR",
:meta_se_bmi_egfr => "META_SE_BMI_eGFR",
:meta_p_bmi_egfr_qval => "META_Q_BMI_eGFR",
:jhs_rho => "JHS_Rho",
:jhs_p => "JHS_P",
:jhs_n => "JHS_N",
:jhs_p_qval => "JHS_Q",
:jhs_rho_bmi => "JHS_Rho_BMI",
:jhs_p_bmi => "JHS_P_BMI",
:jhs_n_bmi => "JHS_N_BMI",
:jhs_p_bmi_qval => "JHS_Q_BMI",
:jhs_rho_bmi_egfr => "JHS_Rho_BMI_eGFR",
:jhs_p_bmi_egfr => "JHS_P_BMI_eGFR",
:jhs_n_bmi_egfr => "JHS_N_BMI_eGFR",
:jhs_p_bmi_egfr_qval => "JHS_Q_BMI_eGFR",
:mesa_rho => "MESA_Rho",
:mesa_p => "MESA_P",
:mesa_n => "MESA_N",
:mesa_p_qval => "MESA_Q",
:mesa_rho_bmi => "MESA_Rho_BMI",
:mesa_p_bmi => "MESA_P_BMI",
:mesa_n_bmi => "MESA_N_BMI",
:mesa_p_bmi_qval => "MESA_Q_BMI",
:mesa_rho_bmi_egfr => "MESA_Rho_BMI_eGFR",
:mesa_p_bmi_egfr => "MESA_P_BMI_eGFR",
:mesa_n_bmi_egfr => "MESA_N_BMI_eGFR",
:mesa_p_bmi_egfr_qval => "MESA_Q_BMI_eGFR",
:her_rho => "HERITAGE_Rho",
:her_p => "HERITAGE_P",
:her_n => "HERITAGE_N",
:her_p_qval => "HERITAGE_Q",
:her_rho_bmi => "HERITAGE_Rho_BMI",
:her_p_bmi => "HERITAGE_P_BMI",
:her_n_bmi => "HERITAGE_N_BMI",
:her_p_bmi_qval => "HERITAGE_Q_BMI")

s2 = rename(s2,s2_names_dict)

sort!(s2,[:Somalogic_ID,:Metabolite_Name])

CSV.write(project_path*tab_out_sub_path*"s2_correlations.csv",s2)


### Supplementary Table 3 - enrichment results
enrichment_joined_file = read_conf(conf,"enrichment_joined_file")
enrichment_joined = CSV.read(project_path*results_sub_path*enrichment_joined_file,DataFrame)
enrichment_joined = rename(enrichment_joined,Dict(:var1 => "SOMA_ID"))
enrichment_joined = leftjoin(enrichment_joined,select(kegg_annot,[:SOMA_ID,:Target_Gene,:Target_Name]), on = :SOMA_ID)

enrichment_joined = select(enrichment_joined,[:SOMA_ID,:Target_Gene,:Target_Name,:class,:es_comb_q])

enrichment_joined.class = replace.(enrichment_joined.class,"var2_" => "")
enrichment_joined.class = replace.(enrichment_joined.class,"Combined_" => "")
enrichment_joined.class = replace.(enrichment_joined.class,"_" => " ")
enrichment_joined.class = lowercase.(enrichment_joined.class)
enrichment_joined.class = replace.(enrichment_joined.class,"organic acids" => "amino & organic acids")

#for i in 1:size(enrichment_joined)[1]
#    t_class = enrichment_joined.class[i][6:end]
#    enrichment_joined.class[i] = t_class
#end

s3 = enrichment_joined

s3 = rename(s3, Dict(:class => "Metabolite_Class", :es_comb_q => "Enrichment_Q"))

CSV.write(project_path*tab_out_sub_path*"s3_enrichment.csv",s3)


### Supplementary Table 7
# analagous to correlation table (by study, all mr methods (no adjustments))

mr_file = read_conf(conf,"mr_file")
mr = CSV.read(project_path*results_sub_path*mr_file,DataFrame)

mr = rename(mr,Dict(Symbol("p.soma_id") => "SOMA_ID", Symbol("m.id") => "ogm_id"))

mr = leftjoin(mr,select(kegg_annot,[:SOMA_ID,:Target_Gene,:Target_Name]), on = :SOMA_ID)
mr = leftjoin(mr,select(ogm_meta,[:ogm_id,:Metabolite_Name]), on = :ogm_id)

s7_names_dict = Dict(
:SOMA_ID => "Somalogic_ID",
Symbol("META_B_ivw") => "META_B_ivw",
Symbol("META_SE_ivw") => "META_SE_ivw",
Symbol("META_P_ivw") => "META_P_ivw",
Symbol("META_ivw_qvalue") => "META_Q_ivw",
Symbol("META_B_maxlik") => "META_B_LIML",
Symbol("META_SE_maxlik") => "META_SE_LIML",
Symbol("META_P_maxlik") => "META_P_LIML",
Symbol("META_maxlik_qvalue") => "META_Q_LIML",
Symbol("META_B_med") => "META_B_median",
Symbol("META_SE_med") => "META_SE_median",
Symbol("META_P_med") => "META_P_median",
Symbol("META_B_med_wt") => "META_B_median_weighted",
Symbol("META_SE_med_wt") => "META_SE_median_weighted",
Symbol("META_P_med_wt") => "META_P_median_weighted",
Symbol("META_B_egg") => "META_B_Egger",
Symbol("META_SE_egg") => "META_SE_Egger",
Symbol("META_P_egg") => "META_P_Egger",
Symbol("META_B_egg_int") => "META_B_Egger_intercept",
Symbol("META_SE_egg_int") => "META_SE_Egger_intercept",
Symbol("META_P_egg_int") => "META_P_Egger_intercept",
Symbol("r.JHS_B_ivw") => "JHS_B_ivw",
Symbol("r.JHS_SE_ivw") => "JHS_SE_ivw",
Symbol("r.JHS_P_ivw") => "JHS_P_ivw",
Symbol("r.JHS_B_maxlik") => "JHS_B_LIML",
Symbol("r.JHS_SE_maxlik") => "JHS_SE_LIML",
Symbol("r.JHS_P_maxlik") => "JHS_P_LIML",
Symbol("r.JHS_B_med") => "JHS_B_median",
Symbol("r.JHS_SE_med") => "JHS_SE_median",
Symbol("r.JHS_P_med") => "JHS_P_median",
Symbol("r.JHS_B_med_wt") => "JHS_B_median_weighted",
Symbol("r.JHS_SE_med_wt") => "JHS_SE_median_weighted",
Symbol("r.JHS_P_med_wt") => "JHS_P_median_weighted",
Symbol("r.JHS_B_egg") => "JHS_B_Egger",
Symbol("r.JHS_SE_egg") => "JHS_SE_Egger",
Symbol("r.JHS_P_egg") => "JHS_P_Egger",
Symbol("r.JHS_B_egg_int") => "JHS_B_Egger_intercept",
Symbol("r.JHS_SE_egg_int") => "JHS_SE_Egger_intercept",
Symbol("r.JHS_P_egg_int") => "JHS_P_Egger_intercept",
Symbol("r.MESA_B_ivw") => "MESA_B_ivw",
Symbol("r.MESA_SE_ivw") => "MESA_SE_ivw",
Symbol("r.MESA_P_ivw") => "MESA_P_ivw",
Symbol("r.MESA_B_maxlik") => "MESA_B_LIML",
Symbol("r.MESA_SE_maxlik") => "MESA_SE_LIML",
Symbol("r.MESA_P_maxlik") => "MESA_P_LIML",
Symbol("r.MESA_B_med") => "MESA_B_median",
Symbol("r.MESA_SE_med") => "MESA_SE_median",
Symbol("r.MESA_P_med") => "MESA_P_median",
Symbol("r.MESA_B_med_wt") => "MESA_B_median_weighted",
Symbol("r.MESA_SE_med_wt") => "MESA_SE_median_weighted",
Symbol("r.MESA_P_med_wt") => "MESA_P_median_weighted",
Symbol("r.MESA_B_egg") => "MESA_B_Egger",
Symbol("r.MESA_SE_egg") => "MESA_SE_Egger",
Symbol("r.MESA_P_egg") => "MESA_P_Egger",
Symbol("r.MESA_B_egg_int") => "MESA_B_Egg_intercept",
Symbol("r.MESA_SE_egg_int") => "MESA_SE_Egg_intercept",
Symbol("r.MESA_P_egg_int") => "MESA_P_Egger_intercept",
Symbol("r.HERITAGE_B_ivw") => "HERITAGE_B_ivw",
Symbol("r.HERITAGE_SE_ivw") => "HERITAGE_SE_ivw",
Symbol("r.HERITAGE_P_ivw") => "HERITAGE_P_ivw",
Symbol("r.HERITAGE_B_maxlik") => "HERITAGE_B_LIML",
Symbol("r.HERITAGE_SE_maxlik") => "HERITAGE_SE_LIML",
Symbol("r.HERITAGE_P_maxlik") => "HERITAGE_P_LIML",
Symbol("r.HERITAGE_B_med") => "HERITAGE_B_median",
Symbol("r.HERITAGE_SE_med") => "HERITAGE_SE_median",
Symbol("r.HERITAGE_P_med") => "HERITAGE_P_median",
Symbol("r.HERITAGE_B_med_wt") => "HERITAGE_B_median_weighted",
Symbol("r.HERITAGE_SE_med_wt") => "HERITAGE_SE_median_weighted",
Symbol("r.HERITAGE_P_med_wt") => "HERITAGE_P_median_weighted",
Symbol("r.HERITAGE_B_egg") => "HERITAGE_B_Egger",
Symbol("r.HERITAGE_SE_egg") => "HERITAGE_SE_Egger",
Symbol("r.HERITAGE_P_egg") => "HERITAGE_P_Egger",
Symbol("r.HERITAGE_B_egg_int") => "HERITAGE_B_Egger_intercept",
Symbol("r.HERITAGE_SE_egg_int") => "HERITAGE_SE_Egger_intercept",
Symbol("r.HERITAGE_P_egg_int") => "HERITAGE_P_Egger_intercept")

mr = rename(mr,s7_names_dict)

s7_cols = [:Somalogic_ID,:Target_Gene,:Target_Name,:Metabolite_Name,:META_B_ivw,:META_SE_ivw,:META_P_ivw,:META_Q_ivw,:META_B_LIML,:META_SE_LIML,:META_P_LIML,:META_Q_LIML,:META_B_median,:META_SE_median,:META_P_median,:META_B_median_weighted,:META_SE_median_weighted,:META_P_median_weighted,:META_B_Egger,:META_SE_Egger,:META_P_Egger,:META_B_Egger_intercept,:META_SE_Egger_intercept,:META_P_Egger_intercept,:JHS_B_ivw,:JHS_SE_ivw,:JHS_P_ivw,:JHS_B_LIML,:JHS_SE_LIML,:JHS_P_LIML,:JHS_B_median,:JHS_SE_median,:JHS_P_median,:JHS_B_median_weighted,:JHS_SE_median_weighted,:JHS_P_median_weighted,:JHS_B_Egger,:JHS_SE_Egger,:JHS_P_Egger,:JHS_B_Egger_intercept,:JHS_SE_Egger_intercept,:JHS_P_Egger_intercept,:MESA_B_ivw,:MESA_SE_ivw,:MESA_P_ivw,:MESA_B_LIML,:MESA_SE_LIML,:MESA_P_LIML,:MESA_B_median,:MESA_SE_median,:MESA_P_median,:MESA_B_median_weighted,:MESA_SE_median_weighted,:MESA_P_median_weighted,:MESA_B_Egger,:MESA_SE_Egger,:MESA_P_Egger,:MESA_B_Egg_intercept,:MESA_SE_Egg_intercept,:MESA_P_Egger_intercept,:HERITAGE_B_ivw,:HERITAGE_SE_ivw,:HERITAGE_P_ivw,:HERITAGE_B_LIML,:HERITAGE_SE_LIML,:HERITAGE_P_LIML,:HERITAGE_B_median,:HERITAGE_SE_median,:HERITAGE_P_median,:HERITAGE_B_median_weighted,:HERITAGE_SE_median_weighted,:HERITAGE_P_median_weighted,:HERITAGE_B_Egger,:HERITAGE_SE_Egger,:HERITAGE_P_Egger,:HERITAGE_B_Egger_intercept,:HERITAGE_SE_Egger_intercept,:HERITAGE_P_Egger_intercept]

s7 = select(mr,s7_cols)

sort!(s7,[:Somalogic_ID,:Metabolite_Name])

CSV.write(project_path*tab_out_sub_path*"s7_mr.csv",s7)



### Supplementary Table 4/5/6 - MR Instruments
# SNP VS METABOLITE GWAS RESULTS
# use only SNPS we have for the study
# include N for exposure

jhs_instr_file = read_conf(conf,"jhs_instr_file")
mesa_instr_file = read_conf(conf,"mesa_instr_file")
her_instr_file = read_conf(conf,"her_instr_file")

instr_names_dict = Dict(
    :A1 => "EFFECT_ALLELE_EXPOSURE",
    :A2 => "OTHER_ALLELE_EXPOSURE",
    :N => "N_EXPOSURE",
    :AF1 => "EAF_EXPOSURE",
    :BETA => "BETA_EXPOSURE",
    :SE => "SE_EXPOSURE",
    :P => "P_EXPOSURE",
    :A1_1 => "EFFECT_ALLELE_OUTCOME",
    :A2_1 => "OTHER_ALLELE_OUTCOME",
    :N_1 => "N_OUTCOME",
    :AF1_1 => "EAF_OUTCOME",
    :BETA_1 => "BETA_OUTCOME",
    :SE_1 => "SE_OUTCOME",
    :P_1 => "P_OUTCOME",
    :exposure => "EXPOSURE",
    :outcome => "OUTCOME")

jhs_instr = CSV.read(project_path*results_sub_path*jhs_instr_file, DataFrame)
filter!(row -> !startswith(row.outcome,"SL"), jhs_instr)
jhs_instr = rename(jhs_instr,instr_names_dict)
jhs_instr = rename(jhs_instr,Dict(:EXPOSURE => "SOMA_ID", :OUTCOME => "jhs_id"))
jhs_meta = filter(row -> !ismissing(row.jhs_id),ogm_meta)[:,[:ogm_id,:refmet_name,:jhs_id]]
jhs_instr = leftjoin(jhs_instr,jhs_meta,on = :jhs_id)
jhs_instr = leftjoin(jhs_instr,select(kegg_annot,[:SOMA_ID,:Target_Gene,:Target_Name]), on = :SOMA_ID)

mr_jhs = filter(row -> !ismissing(row.JHS_B_ivw),mr)
mr_jhs_pairs_set = Set(mr_jhs.Somalogic_ID .* mr_jhs.ogm_id)

filter!(row -> in(row.SOMA_ID * row.ogm_id, mr_jhs_pairs_set),jhs_instr)
jhs_instr = rename(jhs_instr,Dict(:SOMA_ID => "Somalogic_ID_EXPOSURE", :Target_Gene => "Target_Gene_EXPOSURE", :Target_Name => "Target_Name_EXPOSURE", :refmet_name => "Metabolite_Name_OUTCOME"))

s4 = select(jhs_instr,[:SNP,:Somalogic_ID_EXPOSURE,:Target_Gene_EXPOSURE,:Target_Name_EXPOSURE,:EFFECT_ALLELE_EXPOSURE,:OTHER_ALLELE_EXPOSURE,:N_EXPOSURE,:EAF_EXPOSURE,:BETA_EXPOSURE,:SE_EXPOSURE,:P_EXPOSURE,:Metabolite_Name_OUTCOME,:EFFECT_ALLELE_OUTCOME,:OTHER_ALLELE_OUTCOME,:EAF_OUTCOME,:BETA_OUTCOME,:SE_OUTCOME,:P_OUTCOME])

CSV.write(project_path*tab_out_sub_path*"s4_jhs_instruments.csv",s4)



mesa_instr = CSV.read(project_path*results_sub_path*mesa_instr_file, DataFrame)
filter!(row -> !startswith(row.outcome,"SL"), mesa_instr)
mesa_instr = rename(mesa_instr,instr_names_dict)
mesa_instr = rename(mesa_instr,Dict(:EXPOSURE => "SOMA_ID", :OUTCOME => "mesa_id"))
mesa_meta = filter(row -> !ismissing(row.mesa_id),ogm_meta)[:,[:ogm_id,:refmet_name,:mesa_id]]
mesa_instr = leftjoin(mesa_instr,mesa_meta,on = :mesa_id)
mesa_instr = leftjoin(mesa_instr,select(kegg_annot,[:SOMA_ID,:Target_Gene,:Target_Name]), on = :SOMA_ID)

mr_mesa = filter(row -> !ismissing(row.MESA_B_ivw),mr)
mr_mesa_pairs_set = Set(mr_mesa.Somalogic_ID .* mr_mesa.ogm_id)

filter!(row -> in(row.SOMA_ID * row.ogm_id, mr_mesa_pairs_set),mesa_instr)
mesa_instr = rename(mesa_instr,Dict(:SOMA_ID => "Somalogic_ID_EXPOSURE", :Target_Gene => "Target_Gene_EXPOSURE", :Target_Name => "Target_Name_EXPOSURE", :refmet_name => "Metabolite_Name_OUTCOME"))

s5 = select(mesa_instr,[:SNP,:Somalogic_ID_EXPOSURE,:Target_Gene_EXPOSURE,:Target_Name_EXPOSURE,:EFFECT_ALLELE_EXPOSURE,:OTHER_ALLELE_EXPOSURE,:N_EXPOSURE,:EAF_EXPOSURE,:BETA_EXPOSURE,:SE_EXPOSURE,:P_EXPOSURE,:Metabolite_Name_OUTCOME,:EFFECT_ALLELE_OUTCOME,:OTHER_ALLELE_OUTCOME,:EAF_OUTCOME,:BETA_OUTCOME,:SE_OUTCOME,:P_OUTCOME])

CSV.write(project_path*tab_out_sub_path*"s5_mesa_instruments.csv",s5)




her_instr = CSV.read(project_path*results_sub_path*her_instr_file, DataFrame)
filter!(row -> !startswith(row.outcome,"SL"), her_instr)
filter!(row -> !startswith(row.outcome,"QI"), her_instr)
filter!(row -> !startswith(row.outcome,"TF"), her_instr)
her_instr = rename(her_instr,instr_names_dict)
her_instr = rename(her_instr,Dict(:EXPOSURE => "SOMA_ID", :OUTCOME => "heritage_id"))
her_meta = filter(row -> !ismissing(row.heritage_id),ogm_meta)[:,[:ogm_id,:refmet_name,:heritage_id]]
her_instr = leftjoin(her_instr,her_meta,on = :heritage_id)
her_instr = leftjoin(her_instr,select(kegg_annot,[:SOMA_ID,:Target_Gene,:Target_Name]), on = :SOMA_ID)

mr_her = filter(row -> !ismissing(row.HERITAGE_B_ivw),mr)
mr_her_pairs_set = Set(mr_her.Somalogic_ID .* mr_her.ogm_id)

filter!(row -> in(row.SOMA_ID * row.ogm_id, mr_her_pairs_set),her_instr)
her_instr = rename(her_instr,Dict(:SOMA_ID => "Somalogic_ID_EXPOSURE", :Target_Gene => "Target_Gene_EXPOSURE", :Target_Name => "Target_Name_EXPOSURE", :refmet_name => "Metabolite_Name_OUTCOME"))

s6 = select(her_instr,[:SNP,:Somalogic_ID_EXPOSURE,:Target_Gene_EXPOSURE,:Target_Name_EXPOSURE,:EFFECT_ALLELE_EXPOSURE,:OTHER_ALLELE_EXPOSURE,:N_EXPOSURE,:EAF_EXPOSURE,:BETA_EXPOSURE,:SE_EXPOSURE,:P_EXPOSURE,:Metabolite_Name_OUTCOME,:EFFECT_ALLELE_OUTCOME,:OTHER_ALLELE_OUTCOME,:EAF_OUTCOME,:BETA_OUTCOME,:SE_OUTCOME,:P_OUTCOME])

CSV.write(project_path*tab_out_sub_path*"s6_heritage_instruments.csv",s6)