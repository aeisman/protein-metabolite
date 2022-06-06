### compute heritability from master SNP DataFrame

using DataFrames,CSV

try
    include("../analysis/omics_tools.jl")
    include("../analysis/enrichment_tools.jl")
catch
    include("./analysis/omics_tools.jl")
    include("./analysis/enrichment_tools.jl")
end

conf = "./analysis/analysis.conf";
project_path = read_conf(conf,"project_path");
results_sub_path = read_conf(conf,"results_sub_path");
master_snp_file = read_conf(conf,"master_snp_file")

master_snp_df = CSV.read(project_path*results_sub_path*master_snp_file,DataFrame)

function heritability(beta,af1)
    #BETA^2×(2×AF1× (1-AF1)/VAR) where BETA was the beta estimate for the effect allele, AF1 was the allele frequency of the effect allele, and VAR was the variance of the protein residual used for WGS analysis.

    #assume var = 1 (log normalized data)

    var = 1

    heritability = beta^2 * (2 * af1 * (1-af1) / var)

    return heritability
end

function add_heritability_column!(df,beta_col,af1_col,new_col_name)
    if in(new_col_name,Set(names(df)))
        @warn("$(new_col_name) already in df!")
    else
        df[!,Symbol(new_col_name)] .= 0.0
        for i in 1:size(df)[1]
            t_af1 = df[i,Symbol(af1_col)]
            t_beta = df[i,Symbol(beta_col)]
            if t_af1 > 0.0 && !ismissing(t_beta)
                df[i,Symbol(new_col_name)] = heritability(t_beta,t_af1)
            end
        end
    end

    return df
end

add_heritability_column!(master_snp_df,"jhs_beta","jhs_maf","jhs_heritability")
add_heritability_column!(master_snp_df,"mesa_beta","mesa_maf","mesa_heritability")
add_heritability_column!(master_snp_df,"her_beta","her_maf","her_heritability")

function calc_sum_heritability_dict(df,prot_col,select_col,heritability_col)

    prot_set = Set(df[:,Symbol(prot_col)])
    heritability_dict = Dict()
    for prot in prot_set
        prot_df = filter(row -> row[Symbol(prot_col)] == prot && row[Symbol(select_col)] == 1, df)
        heritability_dict[prot] = sum(prot_df[:,Symbol(heritability_col)])
    end

    return heritability_dict
end

jhs_heritability_dict = calc_sum_heritability_dict(master_snp_df,"prot","jhs_select","jhs_heritability")
mesa_heritability_dict = calc_sum_heritability_dict(master_snp_df,"prot","mesa_select","mesa_heritability")
her_heritability_dict = calc_sum_heritability_dict(master_snp_df,"prot","her_select","her_heritability")

function calc_prots_loss(heritability_dict,thresh)
    tot_prots = sum(collect(values(heritability_dict)) .> 0)
    kept_prots = sum(collect(values(heritability_dict)) .> thresh)
    lost_prots = tot_prots - kept_prots
    kept_pct = kept_prots / tot_prots

    return lost_prots, kept_pct
end

thresh_arr = [ 0.05, 0.025, 0.01]
study_dict_arr = [jhs_heritability_dict,mesa_heritability_dict,her_heritability_dict]

for thresh in thresh_arr
    println(thresh)
    for study_dict in study_dict_arr
        (lost_prots,keep_pct) = calc_prots_loss(study_dict,thresh)
        println("$(lost_prots) $(keep_pct)")
    end
    println()
end

function prot_heritability(prot,study_dict_arr)
    println(prot)
    for study_dict in study_dict_arr
        println(study_dict[prot])
    end
end

prots_arr = ["SL000276","SL000668","SL005574"]

for prot in prots_arr
    println()
    prot_heritability(prot,study_dict_arr)
end


jhs_heritability_arr = filtercollect(values(jhs_heritability_dict))


mr = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/mr_filt_sig_final_20220202.csv",DataFrame)
mr_sig = filter(row -> row.META_sig == 3,mr)
mr_sig_prots = Set(mr_sig[:,Symbol("p.soma_id")])

save_prots = Set()
thresh = 0.025
for prot in keys(jhs_heritability_dict)
    if jhs_heritability_dict[prot] > thresh
        push!(save_prots,prot)
    end
end
for prot in keys(mesa_heritability_dict)
    if mesa_heritability_dict[prot] > thresh
        push!(save_prots,prot)
    end
end
for prot in keys(her_heritability_dict)
    if her_heritability_dict[prot] > thresh
        push!(save_prots,prot)
    end
end
intersect(save_prots,mr_sig_prots)


size(filter(row -> in(row[Symbol("p.soma_id")],save_prots),mr_sig))[1]