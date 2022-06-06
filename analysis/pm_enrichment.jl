using CSV, DataFrames, RCall, Distributions, Plots, StatsBase
try
    include("../analysis/omics_tools.jl")
    include("../analysis/enrichment_tools.jl")
catch
    include("./analysis/omics_tools.jl")
    include("./analysis/enrichment_tools.jl")
end

function call_ogea(df,var1,var2,var1_cat,var2_cat,var_p,p,prot_anno,nom,sortby,assume_normal,var2_distr)
    # ogea - omics graph enrichment analysis
    # enrichment analysis modeled after GSEA (PNAS - Subramanian 2005)

    # create protein annotation dict
    prot_anno_dict = Dict()
    for i in 1:size(prot_anno)[1]
        #prot_anno_dict[prot_anno.soma_id[i]] = prot_anno.target_full_name[i]
        prot_anno_dict[prot_anno.soma_id[i]] = split(split(prot_anno.entrez_gene_symbol[i]," ")[1],",")[1]
    end

    t_df = filter(row -> !ismissing(row[Symbol(var_p)]), df) ### remove rows with missing values in var_p column
    t_df = filter(row -> row[Symbol(var1_cat)] == 1 && abs(row[Symbol(var_p)]) <= 1 , t_df)
    sort!(t_df, [Symbol(var_p)])

    if sortby == "abs"
        sort!(t_df, [Symbol("meta_rho")],by = abs, rev = true)
        #sort!(t_df, [Symbol(var_p)],by = abs)
    else
        sort!(t_df, [Symbol("meta_rho")], rev = true)
        #sort!(t_df, [Symbol(var_p)])
    end
    #do ogea for each var1 value
    var1_set = Set(t_df[:,Symbol(var1)])

    std(t_df.meta_rho)

    var_arr = []
    es_arr = []
    max_length = maximum(values(countmap(t_df.var1)))
    ogea_graph_df = DataFrame()

    for var in var1_set
        #var =  "SL003650"
        #var = "SL005574"
        #println(var)
        if haskey(prot_anno_dict,var)
            arr_for_ogea = Matrix(filter(row -> row[Symbol(var1)] == var,t_df)[:,[Symbol(var2_cat),Symbol(var_cor)]])
            var_cum_sum = ogea(arr_for_ogea,p)
            if length(var_cum_sum) < max_length
                while max_length > length(var_cum_sum)
                    push!(var_cum_sum,0)
                end
            end

            ogea_graph_df[:,Symbol(var)] = var_cum_sum
            es = ogea_es(var_cum_sum)
            #println("$(var) $(es)")
            var_arr = vcat(var_arr,var)
            es_arr = vcat(es_arr,es)
        end
    end

    # assemble results into annotated table

    ogea_results = DataFrame(var1 = var_arr, es = es_arr)
    sort!(ogea_results,[:es], rev = true)

    ogea_results[:,:var1_name] .= ""

    for i in 1:size(ogea_results)[1]
        if haskey(prot_anno_dict,ogea_results.var1[i])
            ogea_results.var1_name[i] = prot_anno_dict[ogea_results.var1[i]]
        end
    end

    # determine significance threshold
    n_pos = convert(Int64,round(sum(t_df[:,Symbol(var2_cat)]) / length(var1_set)))
    n_neg = convert(Int64,round(size(t_df)[1] / length(var1_set) - n_pos))

    mu = mean(t_df[:,Symbol(var_cor)])
    sd = std(t_df[:,Symbol(var_cor)])

    cor_data = DataFrame(var_p = t_df[:,Symbol(var_p)],var_cor = t_df[:,Symbol(var_cor)],var2_cat = t_df[:,Symbol(var2_cat)])

    #convert nom to fwes
    #nom = nom / length(var1_set)
    
    #(sig_h,sig_l,es_arr) = ogea_sig_thresh(n_pos, n_neg, 10000, nom, cor_data, 0, 1, sortby)
    # (sig_h,sig_l,es_arr) = ogea_sig_thresh(100,265,100000,0.05,cor_data,1,1,"abs",0)
    # (sig_h,sig_l,es_arr2) = ogea_sig_thresh(163,202,100000,0.05,cor_data,1,1,"abs",0)
    (sig_h,sig_l,es_arr) = ogea_sig_thresh(n_pos, n_neg, 100000, nom, cor_data, p, assume_normal, sortby, var2_distr)
    sort!(es_arr,rev = true)
    
    ogea_results[:,Symbol("sig_greater")] = ogea_results.es .> sig_h[1]

    ogea_results[!,:es_p] .= -999.99
    for i in 1:size(ogea_results)[1]
        #println(i)
        t_es_p = calc_p_from_dist(es_arr,ogea_results.es[i])
        ogea_results.es_p[i] = t_es_p
    end
    #==
    es_norm[!,:es_norm_p] .= -999.99
    es_arr_norm = es_arr / mean(es_arr)
    sort!(es_arr_norm,rev = true)

    es_norm[!,:es_norm_p2] .= -999.99
    es_arr_norm2 = es_arr2 / mean(es_arr2)
    sort!(es_arr_norm2,rev = true)

    for i in 1:length(es_norm.es_norm)
        t_p = calc_p_from_dist(es_arr_norm,es_norm.es_norm[i])
        es_norm.es_norm_p[i] = t_p

        t_p2 = calc_p_from_dist(es_arr_norm2,es_norm.es_norm[i])
        es_norm.es_norm_p2[i] = t_p2
    end
    
    @rput(es_norm)
        R"""
        library(qvalue)
        fdr_level = 0.05
        qvalues = qvalue(es_norm$es_norm_p,fdr.level = fdr_level)$qvalues
        """
    @rget(qvalues)
    sum(qvalues .< 0.05)
    ==#

    @rput(ogea_results)
        R"""
        library(qvalue)
        fdr_level = 0.05
        qvalues = qvalue(ogea_results$es_p,fdr.level = fdr_level, pi0 = 1)$qvalues
        """
    @rget(qvalues)

    ogea_results[!,:es_q] = qvalues

    return(ogea_results, ogea_graph_df, sig_h)
end

#function run_pm_enrichment()

    conf = "./analysis/analysis.conf";
    project_path = read_conf(conf,"project_path");
    results_sub_path = read_conf(conf,"results_sub_path");
    annot_sub_path = read_conf(conf,"annot_sub_path");

    cor_df = CSV.read(project_path*results_sub_path*"all_cor_filt_meta_qval.csv",DataFrame);
    cor_df.var1 = cor_df.PROT;
    cor_df.var2 = cor_df.OGM;

    #var1_remove_set = Set(["SL000277","SL004668","SL004669","SL000424"])
    #var1_remove_set = Set(["SL000277","SL004668","SL004669"])
    #filter!(row -> !in(row.var1,var1_remove_set),cor_df)

    ##  LOAD ANNOTATION INFORMATION

    # PROTEIN ANNOTATIONS
    prot_anno = CSV.read(project_path*annot_sub_path*"kegg_annotated_proteins.csv",DataFrame);

    p_enzyme = CSV.read(project_path*annot_sub_path*"pm_enrichment/enzymes.csv",DataFrame);
    p_plasma_binding = CSV.read(project_path*annot_sub_path*"pm_enrichment/plasma_binding.csv",DataFrame);
    p_transporter_receptor = CSV.read(project_path*annot_sub_path*"pm_enrichment/transporters_receptors.csv",DataFrame);
    p_gpcr = CSV.read(project_path*annot_sub_path*"pm_enrichment/GPCRs.csv",DataFrame);

    p_enzyme_set = Set(filter(row -> startswith(row.soma_id, "SL"), p_enzyme).soma_id);
    p_plasma_binding_set = Set(filter(row -> startswith(row.soma_id, "SL"), p_plasma_binding).soma_id);
    p_transporter_receptor_set = Set(filter(row -> startswith(row.soma_id, "SL"), p_transporter_receptor).soma_id);
    p_gpcr_set = Set(filter(row -> startswith(row.soma_id, "SL"), p_gpcr).soma_id);

    p_transferase_set = Set(filter(row -> !ismissing(row.kegg_enzymes) && startswith(row.kegg_enzymes,"2."),prot_anno).soma_id);
    p_hydrolase_set = Set(filter(row -> !ismissing(row.kegg_enzymes) && startswith(row.kegg_enzymes,"3."),prot_anno).soma_id);

    # add annotations to prot_anno
    prot_anno[!,:p_enzyme] = in.(prot_anno.soma_id,Ref(p_enzyme_set));
    prot_anno[!,:p_plasma_binding] = in.(prot_anno.soma_id,Ref(p_plasma_binding_set));
    prot_anno[!,:p_transporter_receptor] = in.(prot_anno.soma_id,Ref(p_transporter_receptor_set));
    prot_anno[!,:p_gpcr] = in.(prot_anno.soma_id,Ref(p_gpcr_set));

    prot_anno[!,:p_transferase] = in.(prot_anno.soma_id,Ref(p_transferase_set));
    prot_anno[!,:p_hydrolase] = in.(prot_anno.soma_id,Ref(p_hydrolase_set));

    # DO PROTEIN ANNOTATIONS
    cor_df[:,:var1_all] .= 1;
    cor_df[:,:var1_enzyme] = in.(cor_df.var1,Ref(p_enzyme_set));
    cor_df[:,:var1_plasma_binding] = in.(cor_df.var1,Ref(p_plasma_binding_set));
    cor_df[:,:var1_transporter_receptor] = in.(cor_df.var1,Ref(p_transporter_receptor_set));
    cor_df[:,:var1_transferase] = in.(cor_df.var1,Ref(p_transferase_set));
    cor_df[:,:var1_hydrolase] = in.(cor_df.var1,Ref(p_hydrolase_set));
    cor_df[:,:var1_gpcr] = in.(cor_df.var1,Ref(p_gpcr_set));

    cor_df[:,:var1_enzyme_binding_transporter] = in.(cor_df.var1,Ref(union(p_enzyme_set,p_plasma_binding_set,p_transporter_receptor_set)));

    # METABOLITE ANNOTATIONS
    ogm = CSV.read(project_path*annot_sub_path*"ogm_meta_key_anno.csv",DataFrame);
    filter!(row -> in(row.ogm_id,Set(cor_df.OGM)),ogm);

    ogm.pm_class = replace.(ogm.pm_class," " => "_");
    ogm.pm_subclass = replace.(ogm.pm_subclass," " => "_");

    ogm.pm_class = "var2_" .* ogm.pm_class;
    ogm.pm_subclass = "var2_" .* ogm.pm_subclass;

    pm_class_set = Set(ogm.pm_class);
    pm_subclass_set = Set(ogm.pm_subclass);

    for class in pm_class_set
        ogm[!,Symbol(class)] = ogm.pm_class .== class
    end

    for class in pm_subclass_set
        ogm[!,Symbol(class)] = ogm.pm_subclass .== class
    end

    ogm.RefMet_SuperClass = "var2_rm_" .* ogm.RefMet_SuperClass;

    rm_superclass_set = Set(ogm.RefMet_SuperClass)

    for class in rm_superclass_set
        ogm[!,Symbol(class)] = ogm.RefMet_SuperClass .== class
    end

    m_amino_acid_set = Set(filter(row -> row.amino_acid == 1,ogm).ogm_id);
    m_hormone_transmitter_set = Set(filter(row -> row.hormone_transmitter == 1,ogm).ogm_id);
    m_lipid_set = Set(filter(row -> row.lipid == 1,ogm).ogm_id);

    # DO METABOLITE ANNOTATIONS
    for class in union(pm_class_set,pm_subclass_set,rm_superclass_set)
        t_class_set = Set(ogm.ogm_id[Bool.( ogm[:,Symbol(class)] )])
        cor_df[!,Symbol(class)] = in.(cor_df.OGM,Ref(t_class_set))
    end

    cor_df[:,:var2_amino_acid] = in.(cor_df.var2,Ref(m_amino_acid_set));
    cor_df[:,:var2_hormone_transmitter] = in.(cor_df.var2,Ref(m_hormone_transmitter_set));
    cor_df[:,:var2_lipid] = in.(cor_df.var2,Ref(m_lipid_set));



    ## SET OGEA PARAMETERS
    var1 = "var1"
    var2 = "var2"
    var_p = "meta_p"
    var_cor = "meta_rho"
    df = cor_df

    var1_cat_arr = ["var1_all"]
    var2_cat_arr = collect(pm_class_set)
    #var2_cat_arr = ["var2_Combined_Lipids","var2_Organic_acids","var2_Nucleic_acids","var2_Carbohydrates"]
    #var2_cat_arr = collect(union(pm_class_set,rm_superclass_set))
    #var2_cat_arr = ["var2_Organic_acids","var2_rm_Organic acids"]

    p = 1
    nom = 0.05

    #assume_normal_arr = [1,0]
    assume_normal_arr = [1]
    var2_distr_arr = [1,0]
    var2_distr_arr = [1]
    var2_distr_arr = [0]

    sortby = "abs"
    #sortby = "noabs"

    for var2_distr in var2_distr_arr
        for assume_normal in assume_normal_arr
            for var1_cat in var1_cat_arr
                for var2_cat in var2_cat_arr
                    try

                        println("\n\n\ncall_ogea: $(var1),$(var2),$(var1_cat),$(var2_cat),$(var_p),$(p),$(nom),$(sortby),$(assume_normal),$(var2_distr)")
                        t_df = df
                        no_lipids_str = ""
                        if var2_cat != "var2_Combined_Lipids"
                            t_df = filter(row -> row.var2_Combined_Lipids != 1, df)
                            no_lipids_str = "_no_lipids"
                        end
                        println("Run OGEA Analysis...")
                        (ogea_results, ogea_graph_df, sig) = call_ogea(t_df,var1,var2,var1_cat,var2_cat,var_p,p,prot_anno,nom,sortby,assume_normal,var2_distr)
                        println("OGEA Complete!")

                        ## calculate NES
                        mean_es = mean(ogea_results.es)
                        ogea_results[!,:es_norm] = ogea_results.es / mean_es

                        #==
                        ## compute FDR
                        ogea_results[:,:fdr_old] .= false
                        n_sig_fdr_old = size(ogea_results)[1] - convert(Int64,ceil(sum(ogea_results.sig_greater .== 0) / (1/nom - 1) / nom))
                        ogea_results.fdr_old[1:n_sig_fdr_old] .= true

                        ## correct FDR
                        ogea_results[:,:fdr] = ogea_results.es_q .< nom
                        println("FDR complete!")

                        #fdr sig definition
                        #enriched_vars = ogea_results.var1[ogea_results.fdr .== true]
                        #enriched_vars_labels = ogea_results.var1_name[ogea_results.fdr .== true]

                        #fwer sig definition
                        ogea_results[:,:fwer] = ogea_results.es_p .<= nom/length(var2_cat_arr)

                        #sig definition
                        ogea_results[:,:sig] = ogea_results.fwer

                        enriched_vars = ogea_results.var1[ogea_results.sig .== true]
                        enriched_vars_labels = ogea_results.var1_name[ogea_results.sig .== true]

                        ==#

                        println("Writing ogea_results CSV file")
                        ogea_results_fn = project_path*results_sub_path*"enrichment/ogea_results_$(var1)_$(var2)_$(var1_cat)_$(var2_cat)_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal)_var2distr_$(var2_distr)_$(sortby)$(no_lipids_str).csv"
                        CSV.write(ogea_results_fn,ogea_results)
                        println("Writing ogea_graph CSV file")

                        ogea_graph_df_fn = project_path*results_sub_path*"enrichment/ogea_graph_$(var1)_$(var2)_$(var1_cat)_$(var2_cat)_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal)_var2distr_$(var2_distr)_$(sortby)$(no_lipids_str).csv"
                        CSV.write(ogea_graph_df_fn,ogea_graph_df)

                        #==
                        #colors = get_colorRampPalette(color_names[pal],sum(ogea_results.fdr .== true))
                        colors = get_colorRampPalette(color_names[pal],length(enriched_vars))

                        graph_fn = project_path*results_sub_path*"enrichment/ogea_results_$(var1)_$(var2)_$(var1_cat)_$(var2_cat)_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal)_var2distr_$(var2_distr)_$(sortby)$(no_lipids_str).pdf"
                        ogea_graph(ogea_graph_df,enriched_vars,graph_fn,enriched_vars_labels,colors,1)

                        es_graph_fn = project_path*results_sub_path*"enrichment/ogea_results_es_$(var1)_$(var2)_$(var1_cat)_$(var2_cat)_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal)_var2distr_$(var2_distr)_$(sortby)$(no_lipids_str).pdf"
                        ogea_es_graph(ogea_results,es_graph_fn)

                        #pal = pal + 1
                        ==#
                    catch
                        no_lipids_str = ""
                        if var2_cat != "var2_Combined_Lipids"
                            no_lipids_str = "_no_lipids"
                        end
                        println("#####################\n#####################\n#####################\n#####################\n")
                        println("\n call_ogea: $(var1),$(var2),$(var1_cat),$(var2_cat),$(var_p),$(p),$(nom),$(sortby),$(no_lipids_str),$(assume_normal),$(var2_distr) FAILED\n")
                        println("#####################\n#####################\n#####################\n#####################\n")
                    end
                end
            end
        end
    end

    ## combine ogea results
    comb_ogea_results = DataFrame()
    for var2_distr in var2_distr_arr
        for assume_normal in assume_normal_arr
            for var1_cat in var1_cat_arr
                for var2_cat in var2_cat_arr
                    no_lipids_str = ""
                    if var2_cat != "var2_Combined_Lipids"
                        t_df = filter(row -> row.var2_Combined_Lipids != 1, df)
                        no_lipids_str = "_no_lipids"
                    end

                    t_ogea_results_fn = project_path*results_sub_path*"enrichment/ogea_results_$(var1)_$(var2)_$(var1_cat)_$(var2_cat)_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal)_var2distr_$(var2_distr)_$(sortby)$(no_lipids_str).csv"
                    
                    try
                        t_ogea_results = CSV.read(t_ogea_results_fn,DataFrame)
                        t_ogea_results[!,:class] .= var2_cat
                        append!(comb_ogea_results,t_ogea_results)
                    catch
                        @warn("Failed to add $(t_ogea_results_fn)")
                    end

                end
            end
        end
    end

    ## compute comb fdr q values
    @rput(comb_ogea_results)
        R"""
        library(qvalue)
        fdr_level = 0.05
        qvalues = qvalue(comb_ogea_results$es_p,fdr.level = fdr_level, pi0 = 1)$qvalues
        """
    @rget(qvalues)

    comb_ogea_results[!,:es_comb_q] = qvalues
    sum(comb_ogea_results.es_comb_q .< 0.05) ### 882

    comb_ogea_results_sig = filter(row ->row.es_comb_q < nom, comb_ogea_results)

    ## generate enrichment plots of significant results
    #color_names = [["green3","greenyellow"],["dodgerblue4","dodgerblue1"],["darkorchid4","darkorchid1"],["tan4","wheat"],["orange4","orange"],["hotpink4","hotpink"],["darkorchid4","darkorchid1"],["tan4","wheat"],["darkorchid4","darkorchid1"]]
    color_names = [["green3","greenyellow"],["tan4","wheat"],["tan4","wheat"],["tan4","wheat"],["tan4","wheat"],["dodgerblue4","dodgerblue1"],["darkorchid4","darkorchid1"],["hotpink4","hotpink"],["tan4","wheat"],["orange4","orange"],["hotpink4","hotpink"],["darkorchid4","darkorchid1"],["tan4","wheat"],["darkorchid4","darkorchid1"]]
    pal = 1

    for var2_distr in var2_distr_arr
        for assume_normal in assume_normal_arr
            for var1_cat in var1_cat_arr
                for var2_cat in var2_cat_arr
                    no_lipids_str = ""
                    if var2_cat != "var2_Combined_Lipids"
                        t_df = filter(row -> row.var2_Combined_Lipids != 1, df)
                        no_lipids_str = "_no_lipids"
                    end

                    t_ogea_results_fn = project_path*results_sub_path*"enrichment/ogea_results_$(var1)_$(var2)_$(var1_cat)_$(var2_cat)_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal)_var2distr_$(var2_distr)_$(sortby)$(no_lipids_str).csv"
                    t_ogea_graph_df_fn = project_path*results_sub_path*"enrichment/ogea_graph_$(var1)_$(var2)_$(var1_cat)_$(var2_cat)_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal)_var2distr_$(var2_distr)_$(sortby)$(no_lipids_str).csv"
                    try
                        t_ogea_results = CSV.read(t_ogea_results_fn,DataFrame)
                        t_ogea_graph_df = CSV.read(t_ogea_graph_df_fn,DataFrame)

                        println("files loaded")

                        enriched_vars = filter(row -> row.class == var2_cat,comb_ogea_results_sig).var1
                        enriched_vars_labels = filter(row -> row.class == var2_cat,comb_ogea_results_sig).var1_name

                        println("enriched vars defined")

                        colors = get_colorRampPalette(color_names[pal],length(enriched_vars))

                        println("colors defined")

                        graph_fn = project_path*results_sub_path*"enrichment/ogea_results_$(var1)_$(var2)_$(var1_cat)_$(var2_cat)_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal)_var2distr_$(var2_distr)_$(sortby)$(no_lipids_str)_es_comb_q.pdf"
                        ogea_graph(t_ogea_graph_df,enriched_vars,graph_fn,enriched_vars_labels,colors,1)
                        pal = pal + 1
                        println("graph 1 generated")

                        
                        es_graph_fn = project_path*results_sub_path*"enrichment/ogea_results_es_$(var1)_$(var2)_$(var1_cat)_$(var2_cat)_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal)_var2distr_$(var2_distr)_$(sortby)$(no_lipids_str)_es_comb_q.pdf"
                        t_ogea_results[!,:sig] .= false
                        println("defined sig var")
                        
                        for i in 1:size(t_ogea_results)[1]
                            if in(t_ogea_results.var1[i],Set(enriched_vars))
                                t_ogea_results.sig[i] = true
                            end
                        end
                        println("defined sig var")
                        ogea_es_graph(t_ogea_results,es_graph_fn)
                        

                        println("graph 2 generated")
                    catch
                        @warn("Failed to generate $(t_ogea_graph_df_fn) graph")
                    end
                end
            end
        end
    end

    CSV.write(project_path*results_sub_path*"enrichment/comb_enrichment_$(var_p)_$(p)_$(nom)_assume_normal_$(assume_normal_arr[1])_var2distr_$(var2_distr_arr[1])_$(sortby)_no_lipids_es_comb_q.csv",comb_ogea_results)



    # plot distributions of metabolite rho'S

    for var2_cat in var2_cat_arr

        cat_rho = filter(row -> row[Symbol(var2_cat)], cor_df).meta_rho
        fh = project_path*results_sub_path*"/enrichment/$(var2_cat)_rho_distribution.pdf"
        @rput(cat_rho, fh, var2_cat)
        R"""
            pdf(fh, width = 10, height = 10)
            #par(mar=c(2, 3, 1, 1),pin = c(6.75,2.5))
            #mean_rho = mean(cat_rho)
            #sd_rho = sd(cat_rho)
            hist(cat_rho,xlab = var2_cat,xlim=c(-1,1))
            dev.off()
        """
        #sleep(5)
    end

#end

#run_pm_enrichment()
#==
Need to look at the numbers of each metabolite class that are included depending on distr1 vs distr0
Look at placement of ACY1

# julia analysis/combine_files.jl "/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/enrichment" csv
# need to remove extra header rows and create the metabolite column

comb = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/enrichment/COMB_ogea_results_var1_var2_var1_all_meta_p_1_0.05_assume_normal_1_var2distr_0_abs_no_lipids.csv",DataFrame)


@rput(comb)
    R"""
    library(qvalue)
    fdr_level = 0.05
    qvalues = qvalue(comb$es_p,fdr.level = fdr_level)$qvalues
    """
@rget(qvalues)

comb[!,:es_all_q] = qvalues
sum(comb.es_all_q .< 0.05)

### combine files and determine global normalized enrichment scores
# normalize to mean of 1 and standard deviation of 0.1
# standard deviation and mean need to be of the simulated distribution
for var2_distr in var2_distr_arr
    for assume_normal in assume_normal_arr
        for var1_cat in var1_cat_arr
            for var2_cat in var2_cat_arr

                t_df = filter(row -> !ismissing(row[Symbol(var_p)]), df) ### remove rows with missing values in var_p column
                t_df = filter(row -> row[Symbol(var1_cat)] == 1 && abs(row[Symbol(var_p)]) <= 1 , t_df)
                sort!(t_df, [Symbol(var_p)])

                if sortby == "abs"
                    sort!(t_df, [Symbol("meta_rho")],by = abs, rev = true)
                    #sort!(t_df, [Symbol(var_p)],by = abs)
                else
                    sort!(t_df, [Symbol("meta_rho")], rev = true)
                    #sort!(t_df, [Symbol(var_p)])
                end
                cor_data = DataFrame(var_p = t_df[:,Symbol(var_p)],var_cor = t_df[:,Symbol(var_cor)],var2_cat = t_df[:,Symbol(var2_cat)])
                (sig_h,sig_l,es_arr) = ogea_sig_thresh(100,265,100000,0.05,cor_data,1,1,"abs",0)

            end
        end
    end
end
==#

#==
Testing updated enrichment script
e_joined = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/enrichment_joined_final_20220408.csv",DataFrame);
lipids = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/enrichment/ogea_results_var1_var2_var1_all_var2_Combined_Lipids_meta_p_1_0.05_abs.csv",DataFrame);
aa = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/enrichment/ogea_results_var1_var2_var1_all_var2_Organic_acids_meta_p_1_0.05_abs.csv",DataFrame);
carbs = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/enrichment/ogea_results_var1_var2_var1_all_var2_Carbohydrates_meta_p_1_0.05_abs.csv",DataFrame);
nitro = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/enrichment/ogea_results_var1_var2_var1_all_var2_Organic_nitrogen_compounds_meta_p_1_0.05_abs.csv",DataFrame);
cyclics = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/enrichment/ogea_results_var1_var2_var1_all_var2_Organoheterocyclic_compounds_meta_p_1_0.05_abs.csv",DataFrame);
nucleic = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/enrichment/ogea_results_var1_var2_var1_all_var2_Nucleic_acids_meta_p_1_0.05_abs.csv",DataFrame);

e_joined_lip_sig = Set(filter(row -> row.Combined_Lipids_fdr == true,e_joined).var1);
e_joined_aa_sig = Set(filter(row -> row.Organic_acids_fdr == true,e_joined).var1);
e_joined_car_sig = Set(filter(row -> row.Carbohydrates_fdr == true,e_joined).var1);
e_joined_nit_sig = Set(filter(row -> row.Organic_nitrogen_compounds_fdr == true,e_joined).var1);
e_joined_cyc_sig = Set(filter(row -> row.Organoheterocyclic_compounds_fdr == true,e_joined).var1);
e_joined_nuc_sig = Set(filter(row -> row.Nucleic_acids_fdr == true,e_joined).var1);

lip_sig = Set(filter(row -> row.fdr == true,lipids).var1);
aa_sig = Set(filter(row -> row.fdr == true,aa).var1);
car_sig = Set(filter(row -> row.fdr == true,carbs).var1);
nit_sig = Set(filter(row -> row.fdr == true,nitro).var1);
cyc_sig = Set(filter(row -> row.fdr == true,cyclics).var1);
nuc_sig = Set(filter(row -> row.fdr == true,nucleic).var1);

length(e_joined_lip_sig)
length(intersect(e_joined_lip_sig,lip_sig))

length(e_joined_aa_sig)
length(intersect(e_joined_aa_sig,aa_sig))

length(e_joined_car_sig)
length(intersect(e_joined_car_sig,car_sig))

length(e_joined_nit_sig)
length(intersect(e_joined_nit_sig,nit_sig))

length(e_joined_cyc_sig)
length(intersect(e_joined_cyc_sig,cyc_sig))

length(e_joined_nuc_sig)
length(intersect(e_joined_nuc_sig,nuc_sig))
==#