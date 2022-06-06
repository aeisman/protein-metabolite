using DataFrames, CSV, StatsBase, RCall
try
    include("../analysis/omics_tools.jl")
catch
    include("./analysis/omics_tools.jl")
end

function clean_cor_df!(cor_df,met_anno)
    filter!(row -> !ismissing(row.meta_p_qval),cor_df) ## should be none
    filter!(row -> !endswith(row.cor_id,r"[AB]"),cor_df) ## remove duplicate A/B soma_ids from heritage
    filter!(row -> !endswith(row.cor_id,"A_1"),cor_df)

    cor_df_mets = Set(cor_df.OGM)
    met_anno_mets = Set(met_anno.ogm_id)

    met_set = intersect(cor_df_mets,met_anno_mets)

    filter!(row -> in(row.OGM,met_set), cor_df) # only keep mets that are in the meta-anno df (reduced for duplicate across platform)
end

function qvalue(p_arr, pi_val)
    @rput(p_arr,pi_val)
    R"""
        require('qvalue')
        p_arr[p_arr == -999.99] = NA

        if (sum(is.na(p_arr) < length(p_arr))){
            q_arr = qvalue(p_arr, pi0 = pi_val)$qvalues
        } else {
            q_arr = NA
        }
    """
    @rget(q_arr)

    return q_arr
end

function combine_cor(studies_arr)
    conf = "./analysis/analysis.conf" 
    project_path = read_conf(conf,"project_path")
    adj_str_arr = ["","_bmi","_bmi_egfr"]

    all_cor_filt = DataFrame()

    studies_arr = ["jhs","mesa","her"]
    for study in studies_arr
        print("Loading $(study) Correlation Data...")
        #study = "her"
        project_sub_out_path = read_conf(conf,study*"_data_sub_out_path");
        study_cor_filt_fh = project_path*project_sub_out_path*study*"_cor_filt.csv"
        study_cor_filt = CSV.read(study_cor_filt_fh,DataFrame)
        
        study_cor_filt[!,:drop] .= 0

        study_exists_dict = Dict()
        for i in 1:size(study_cor_filt)[1]
            cor_id = study_cor_filt.cor_id[i]
            if haskey(study_exists_dict,cor_id)
                study_cor_filt.drop[i] = 1
            else
                study_exists_dict[cor_id] = 1
            end
        end
        filter!(row -> row.drop == 0, study_cor_filt)

        if size(all_cor_filt) == (0,0)
            all_cor_filt = study_cor_filt
        else
            all_cor_filt = outerjoin(all_cor_filt, study_cor_filt, on = :cor_id, makeunique = true)
        end
        println("Done!!!")
    end

    # generate PROT and OGM column from cor_id (would not work for protein - protein)
    all_cor_filt[:,:OGM] .= ""
    for i in 1:size(all_cor_filt)[1]
        if contains(all_cor_filt.cor_id[i],"OGM")
            all_cor_filt.OGM[i] = match(r"OGM[0-9]+",all_cor_filt.cor_id[i]).match
        end
    end
    all_cor_filt[:,:PROT] .= ""
    for i in 1:size(all_cor_filt)[1]
        if contains(all_cor_filt.cor_id[i],"SL")
            all_cor_filt.PROT[i] = match(r"SL[0-9]+",all_cor_filt.cor_id[i]).match
        end
    end
    # remove HCE rows
    filter!(row -> row.PROT != "",all_cor_filt)

    # only keep proteins that are in jhs or mesa and start with SL
    keep_prots = union(Set(all_cor_filt.jhs_var1_alt),Set(all_cor_filt.mesa_var1_alt))
    filter!(row -> in(row.PROT,keep_prots),all_cor_filt)

    # Remove APOE variants ?,"SL000424"
    remove_prots = Set(["SL000277","SL004668","SL004669"])
    filter!(row -> !in(row.PROT,remove_prots), all_cor_filt)

    # Remove METS (control and drugs)  "atenolol_binary", "Atorvastatin", "Valsartan", "metoprolol_binary", "Metronidazole", "metformin_binary", "Metoprolol", "Acetaminophen", "Simvastatin", "Quinine", "valine-d8", "hydroxycotinine_binary", "Salicylic acid", "Verapamil", "metronidazole_binary", "Lisinopril", "Pravastatin", "Oxypurinol", "phenylalanine-d8", "Metformin", "Warfarin"
 
    exclude_set = Set(["OGM313241","OGM317024","OGM317714","OGM316680","OGM317833","OGM318901","OGM311526","OGM316665","OGM315324","OGM311467","OGM317508","OGM318240","OGM311785","OGM312556","OGM316580","OGM314884","OGM316961","OGM310940","OGM312106","OGM319331","OGM317197","OGM316852"])
    remove_mets = Set(exclude_set)
    filter!(row -> !in(row.OGM,remove_mets), all_cor_filt)


    ### Do Correlation META-ANALYSIS

    println("Doing Correlation Meta-Analysis...")

    all_cor_filt_meta = DataFrame()

    total = size(all_cor_filt)[1]
    for i in 1:size(all_cor_filt)[1]
        row = DataFrame(all_cor_filt[i,:])
        #if i > 471740
            if rem(i,round(0.05*total)) == 0
                pct = i/total*100
                pct_str = @sprintf("%.1f",pct)
                println("$i of $total | $pct_str%")
            end

            for adj_str in adj_str_arr
                #println("adj $(adj_str)")
                r = []
                n = []
                labels = []
                for study in studies_arr
                    if study == "her" && adj_str == "_bmi_egfr" ## skip heritage for egfr adjustment
                    else
                        r_study = row[1,Symbol(study*"_rho"*adj_str)]
                        n_study = row[1,Symbol(study*"_n"*adj_str)]
                        if !ismissing(r_study) ## check for missing values (e.g. protein/metabolite not available for particular study)
                            if size(r) == (0,)
                                r = r_study
                                n = n_study
                            else
                                r = vcat(r,r_study)
                                n = vcat(n,n_study)
                            end
                            push!(labels,study)
                        end
                    end
                end

                #println("assemble df for meta")

                df_for_meta = DataFrame(r = r, n = n, labels = labels)
                filter!(row -> !isnan(row.r) ,df_for_meta)
                #println("filtered!")
                if size(df_for_meta)[1] > 0
                    #println("attempt metacor")
                    @rput(df_for_meta)
                    R"""
                        require('metacor')
                        meta_op = metacor.OP(df_for_meta$r,df_for_meta$n,df_for_meta$labels,plot=FALSE)
                    """
                    @rget(meta_op)

                    row[!,Symbol("meta_rho"*adj_str)] .= meta_op[:G_mean]
                    row[!,Symbol("meta_p"*adj_str)] .= meta_op[:p]
                    row[!,Symbol("meta_se"*adj_str)] .= meta_op[:G_mean_se]
                else
                    #println("assign missing values")
                    row[!,Symbol("meta_rho"*adj_str)] .= missing
                    row[!,Symbol("meta_p"*adj_str)] .= missing
                    row[!,Symbol("meta_se"*adj_str)] .= missing
                    #println("missing values assigned")
                end

            end

            #println("Append!")
            allowmissing!(all_cor_filt_meta)
            append!(all_cor_filt_meta,DataFrame(row))
        #end
        
    end
#==

    # create placeholder values
    all_cor_filt[!,:meta_rho] .= -999.99;
    all_cor_filt[!,:meta_p] .= -999.99;
    all_cor_filt[!,:meta_se] .= -999.99;

    all_cor_filt[!,:meta_rho_bmi] .= -999.99;
    all_cor_filt[!,:meta_p_bmi] .= -999.99;
    all_cor_filt[!,:meta_se_bmi] .= -999.99;

    all_cor_filt[!,:meta_rho_bmi_egfr] .= -999.99;
    all_cor_filt[!,:meta_p_bmi_egfr] .= -999.99;
    all_cor_filt[!,:meta_se_bmi_egfr] .= -999.99;

    total = size(all_cor_filt)[1]
    
    for i in 1:size(all_cor_filt)[1]

        if rem(i,round(0.05*total)) == 0
            pct = i/total*100
            pct_str = @sprintf("%.1f",pct)
            println("$i of $total | $pct_str%")
        end

        #println()
        #println(i)
        for adj_str in adj_str_arr
            #println(adj_str)
            r = []
            n = []
            labels = []
            for study in studies_arr
                if study == "her" && adj_str == "_bmi_egfr" ## skip heritage for egfr adjustment
                else
                    r_study = all_cor_filt[i,Symbol(study*"_rho"*adj_str)]
                    n_study = all_cor_filt[i,Symbol(study*"_n"*adj_str)]
                    if !ismissing(r_study) ## check for missing values (e.g. protein/metabolite not available for particular study)
                        if size(r) == (0,)
                            r = r_study
                            n = n_study
                        else
                            r = vcat(r,r_study)
                            n = vcat(n,n_study)
                        end
                        push!(labels,study)
                    end
                end
            end

            df_for_meta = DataFrame(r = r, n = n, labels = labels)
            filter!(row -> !isnan(row.r) ,df_for_meta)

            if size(df_for_meta)[1] > 0
                @rput(df_for_meta)
                R"""
                    require('metacor')
                    meta_op = metacor.OP(df_for_meta$r,df_for_meta$n,df_for_meta$labels,plot=FALSE)
                """
                @rget(meta_op)

                all_cor_filt[i,Symbol("meta_rho"*adj_str)] = meta_op[:G_mean]
                all_cor_filt[i,Symbol("meta_p"*adj_str)] = meta_op[:p]
                all_cor_filt[i,Symbol("meta_se"*adj_str)] = meta_op[:G_mean_se]
            end

        end

    end
==#
    # remove rows without results
    filter!(row -> !ismissing(row.meta_p), all_cor_filt_meta);
    filter!(row -> !isnan(row.meta_p), all_cor_filt_meta);
    filter!(row -> row.meta_p != -999.99, all_cor_filt_meta);

    println("Done!!!")

    ### Calculate by-protein META qvalues

    println("Calculating by-protein qvalues...")

    cor_prots = collect(Set(all_cor_filt_meta.PROT))
    all_cor_filt_meta_qval = DataFrame()
    
    for i in 1:length(cor_prots)
        prot = cor_prots[i]
        println("$(i) $(prot)")
        prot_df = filter(row -> row.PROT == prot, all_cor_filt_meta)

        for adj_str in adj_str_arr
            prot_df[:,Symbol("meta_p"*adj_str*"_qval")] = qvalue(prot_df[:,Symbol("meta_p"*adj_str)],1)
                
            for study in studies_arr
                if study == "her" && adj_str == "_bmi_egfr" ## skip heritage for egfr adjustment
                else
                    if sum(ismissing.(prot_df[:,Symbol(study*"_p"*adj_str)])) < length(prot_df[:,Symbol(study*"_p"*adj_str)])
                        prot_df[:,Symbol(study*"_p"*adj_str*"_qval")] = qvalue(prot_df[:,Symbol(study*"_p"*adj_str)],1)
                    else
                        prot_df[:,Symbol(study*"_p"*adj_str*"_qval")] .= missing
                    end
                end
            end
        end

        append!(all_cor_filt_meta_qval,prot_df)

    end

    println("Done!!!")

    annot_sub_path = read_conf(conf,"annot_sub_path")
    met_anno_fh = project_path*annot_sub_path*"ogm_meta_key_anno.csv"
    met_anno = CSV.read(met_anno_fh, DataFrame)

    clean_cor_df!(all_cor_filt_meta_qval,met_anno)

    results_sub_path = read_conf(conf,"results_sub_path")
    all_cor_filt_meta_qval_fh = project_path*results_sub_path*"all_cor_filt_meta_qval.csv"

    CSV.write(all_cor_filt_meta_qval_fh,all_cor_filt_meta_qval)

end

conf = "./analysis/analysis.conf" 
project_path = read_conf(conf,"project_path")

studies_arr = ["jhs","mesa","her"]
combine_cor(studies_arr)