# Determine which pairwise protein-metabolite relationships to test with MR and filter MR results
# Main input is joined output from neo4j database

using DataFrames,CSV,StatsBase,RCall
try
    include("../analysis/omics_tools.jl")
catch
    include("./analysis/omics_tools.jl")
end

function mr_fdr(mr_df,var_p)
    df_t = rename(mr_df,Symbol(var_p) => :mr_p_var)
    t_q = ones(size(df_t)[1])
    try
        @rput(df_t)
        R"""
            library(qvalue)
            t_q = qvalue(df_t$mr_p_var)$qvalues
        """
        @rget(t_q)
    catch
        println("Using Benhamini-Hoshberg Method, pi0 = 1")
        @rput(df_t)
        R"""
            library(qvalue)
            t_q = qvalue(df_t$mr_p_var, pi0 = 1)$qvalues
        """
        @rget(t_q)
    end

    q = []
    j = 1

    for i in 1:length(df_t.mr_p_var)
        if !ismissing(df_t.mr_p_var[i])
            push!(q,t_q[j])
            j = j + 1
        else
            push!(q,missing)
        end
    end

    return q
end

function call_mr_fdr_by_var(mr_df,var_p,byvar)
    mr_fdr_df = DataFrame()
    byvar_set = Set(mr_df[:,Symbol(byvar)])

    for var in byvar_set
        println(var)
        t_mr_filt = filter(row -> row[Symbol(byvar)] == var, mr_df)

        if size(t_mr_filt)[1] > 1
            qvalues = mr_fdr(t_mr_filt,var_p)
            t_mr_filt[!,:META_sig_qvalue] = qvalues
            #println("assigned qvalues")
            
        elseif size(t_mr_filt)[1] == 1
            @warn("Only one row!!!")
            t_mr_filt[!,:META_sig_qvalue] = t_mr_filt.META_P_ivw

        else
            @warn("Less than 1 row!!!")
        end
    
        mr_fdr_df = vcat(mr_fdr_df,t_mr_filt)

    end

    return mr_fdr_df
end

conf = "./analysis/analysis.conf";
project_path = read_conf(conf,"project_path");
results_sub_path = read_conf(conf,"results_sub_path");

cor_df = CSV.read(project_path*results_sub_path*"all_cor_filt_meta_qval.csv",DataFrame);
mr_df = CSV.read(project_path*results_sub_path*"/mr/mr_pm_meta_20220429.csv",DataFrame);

mr_df[!,:cor_id] .= ""
## add cor_id to mr df
for i in 1:size(mr_df)[1]
    cor_id_arr = sort([mr_df[i,Symbol("p.soma_id")],mr_df[i,Symbol("m.id")]])
    mr_df.cor_id[i] = join(cor_id_arr)
end


filter!(row -> !ismissing(row.meta_p_qval),cor_df)
### remove APOE variant proteins
prot_remove_set = Set(["SL000277","SL004668","SL004669","SL000424"])
filter!(row -> !in(row.PROT,prot_remove_set),cor_df)
# remove proteins that did not have instruments in any study (and therefore don't have MR results)
filter!(row -> !in(row[Symbol("p.soma_id")],prot_remove_set),mr_df)

### Determine pairwise cor_id's to bring forward from correlation analysis
sig = 0.05
cor_df_sig_set = Set(filter(row -> row.meta_p_qval < sig, cor_df).cor_id)

### Determine pairwise cor_id's to bring forward from enrichment analysis
enrichment_dir = project_path*results_sub_path*"enrichment/final_20220421/"

e_joined = CSV.read(enrichment_dir*"comb_enrichment_meta_p_1_0.05_assume_normal_1_var2distr_0_abs_no_lipids_es_comb_q.csv",DataFrame)

# remove "var2_ from class names
e_joined.class = replace.(e_joined.class,"var2_" => "")
filter!(row -> row.es_comb_q < sig, e_joined)

annot_sub_path = read_conf(conf,"annot_sub_path");
ogm = CSV.read(project_path*annot_sub_path*"ogm_meta_key_anno.csv",DataFrame);
ogm.pm_class = replace.(ogm.pm_class," " => "_")
filter!(row -> in(row.ogm_id,Set(cor_df.OGM)),ogm);

sort!(ogm,[:RefMet_SuperClass,:RefMet_Main_Class])

mets = sort(collect(Set(cor_df.OGM)))
filter!(row -> in(row.ogm_id, Set(mets)), ogm)

mets_anno = sort(collect(Set(ogm.ogm_id)))


### Generate set of cor_id's that correspond to enrichment analysis results
#e_joined_vars_to_check = names(e_joined)[endswith.(names(e_joined),"es_q")]
e_joined_sig_set = Set()

for i in 1:size(e_joined)[1]
    t_prot = e_joined.var1[i]
    t_mets = filter(row -> row.pm_class == e_joined.class[i],ogm).ogm_id
    for met in t_mets
        t_cor_id = join(sort([t_prot,met]))
        push!(e_joined_sig_set,t_cor_id)
    end
end

cor_sig_enriched_union = union(cor_df_sig_set,e_joined_sig_set)

cor_sig_enriched_union_df = DataFrame(cor_id = [], cor_sig = [], enriched = [])
for cor_id in cor_sig_enriched_union
    cor_sig = in(cor_id,cor_df_sig_set)
    enriched = in(cor_id,e_joined_sig_set)
    t_df = DataFrame(cor_id = cor_id, cor_sig = cor_sig, enriched = enriched)
    cor_sig_enriched_union_df = vcat(cor_sig_enriched_union_df,t_df)
end

CSV.write(project_path*results_sub_path*"cor_sig_enriched_union.csv", cor_sig_enriched_union_df)

mr_cor_filt = filter(row -> in(row.cor_id,cor_df_sig_set),mr_df);
mr_enrich_filt = filter(row -> in(row.cor_id,e_joined_sig_set),mr_df);
mr_cor_enrich_filt = filter(row -> in(row.cor_id,cor_sig_enriched_union),mr_df);

### compute fdr_sig5/10 by protein for mr_filt
mr_cor_filt_qvalue = call_mr_fdr_by_var(mr_cor_filt, "META_P_ivw", "p.soma_id");
rename!(mr_cor_filt_qvalue,Dict(:META_sig_qvalue => :META_ivw_qvalue));
mr_cor_filt_qvalue = call_mr_fdr_by_var(mr_cor_filt_qvalue, "META_P_maxlik", "p.soma_id");
rename!(mr_cor_filt_qvalue,Dict(:META_sig_qvalue => :META_maxlik_qvalue));

mr_enrich_filt_qvalue = call_mr_fdr_by_var(mr_enrich_filt, "META_P_ivw", "p.soma_id");
rename!(mr_enrich_filt_qvalue,Dict(:META_sig_qvalue => :META_ivw_qvalue));
mr_enrich_filt_qvalue = call_mr_fdr_by_var(mr_enrich_filt_qvalue, "META_P_maxlik", "p.soma_id");
rename!(mr_enrich_filt_qvalue,Dict(:META_sig_qvalue => :META_maxlik_qvalue));

mr_cor_enrich_filt_qvalue = call_mr_fdr_by_var(mr_cor_enrich_filt, "META_P_ivw", "p.soma_id");
rename!(mr_cor_enrich_filt_qvalue,Dict(:META_sig_qvalue => :META_ivw_qvalue));
mr_cor_enrich_filt_qvalue = call_mr_fdr_by_var(mr_cor_enrich_filt_qvalue, "META_P_maxlik", "p.soma_id");
rename!(mr_cor_enrich_filt_qvalue,Dict(:META_sig_qvalue => :META_maxlik_qvalue));


sum(mr_cor_filt_qvalue.META_ivw_qvalue .< sig)
sum(mr_enrich_filt_qvalue.META_ivw_qvalue .< sig)
sum(mr_cor_enrich_filt_qvalue.META_ivw_qvalue .< sig)

mr_cor_filt_qvalue_sig = sort(filter(row -> row.META_ivw_qvalue < sig, mr_cor_filt_qvalue),:META_ivw_qvalue)
mr_enrich_filt_qvalue_sig = sort(filter(row -> row.META_ivw_qvalue < sig, mr_enrich_filt_qvalue),:META_ivw_qvalue)
mr_cor_enrich_filt_qvalue_sig = sort(filter(row -> row.META_ivw_qvalue < sig, mr_cor_enrich_filt_qvalue),:META_ivw_qvalue)

sort(countmap(mr_cor_filt_qvalue_sig[:,Symbol("p.name")]))
sort(countmap(mr_enrich_filt_qvalue_sig[:,Symbol("p.name")]))
sort(countmap(mr_cor_enrich_filt_qvalue_sig[:,Symbol("p.name")]))

### annotate final results with q values from only cor and only enrichment methods of pruning MR results before FDR calculations
cor_only_qvalue = []
enrich_only_qvalue = []

for i in 1:size(mr_cor_enrich_filt_qvalue)[1]
    t_cor_id = mr_cor_enrich_filt_qvalue.cor_id[i]

    t_cor_only = filter(row -> row.cor_id == t_cor_id,mr_cor_filt_qvalue)
    t_enrich_only = filter(row -> row.cor_id == t_cor_id,mr_enrich_filt_qvalue)

    if size(t_cor_only)[1] == 1
        push!(cor_only_qvalue,t_cor_only.META_ivw_qvalue[1])
    else
        push!(cor_only_qvalue,missing)
    end

    if size(t_enrich_only)[1] == 1
        push!(enrich_only_qvalue,t_enrich_only.META_ivw_qvalue[1])
    else
        push!(enrich_only_qvalue,missing)
    end


    if (i % 10000) == 0
        println(i)
    end
end

mr_cor_enrich_filt_qvalue[!,:cor_only_qvalue] = cor_only_qvalue
mr_cor_enrich_filt_qvalue[!,:enrich_only_qvalue] = enrich_only_qvalue

CSV.write(project_path*results_sub_path*"mr_pm_meta_20220429_qval.csv",mr_cor_enrich_filt_qvalue)
