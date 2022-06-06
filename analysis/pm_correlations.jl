using DataFrames,CSV,StatsBase,RCall
try
    include("../analysis/omics_tools.jl")
catch
    include("./analysis/omics_tools.jl")
end

function pm_correlations(studies_arr)
    conf = "./analysis/analysis.conf"
    project_path = read_conf(conf,"project_path")

    for study in studies_arr
        println(study)
        study_suffix = ""
        if study == "her"
            study_suffix = "itage"
        end

        project_sub_out_path = read_conf(conf,study*"_data_sub_out_path");
        project_sub_path = read_conf(conf,study*"_data_sub_path");
        project_sub_path_files_arr = readdir(project_path*project_sub_path);

        prot_names = CSV.read(project_path*project_sub_path*study*"_prot_names.csv", DataFrame, header = false).Column1;
        met_names = CSV.read(project_path*project_sub_path*study*"_met_names.csv", DataFrame, header = false).Column1;

        all_names = vcat(prot_names,met_names)

        ogm_sub_path = read_conf(conf,"ogm_sub_path")
        ogm_meta_fh = project_path*ogm_sub_path*"ogm_meta_key.csv"
        ogm_key_df = CSV.read(ogm_meta_fh,DataFrame);

        filter!(row -> !ismissing(row[Symbol(study*study_suffix*"_platform")]), ogm_key_df)
        ogm_key_ami_df = filter(row -> row[Symbol(study*study_suffix*"_platform")] == "amide", ogm_key_df);
        ogm_key_hil_df = filter(row -> row[Symbol(study*study_suffix*"_platform")] == "hilic", ogm_key_df);

        # only keep met names that are in the ogm_key file
        met_names = collect(intersect(Set(met_names),Set(ogm_key_df[:,Symbol(study*study_suffix*"_name")])))

        study_files_arr = project_sub_path_files_arr[startswith.(project_sub_path_files_arr,study*"_join")];

        ### LOAD 5 to 1.3K Conversion Dictionary
        annot_sub_path = read_conf(conf,"annot_sub_path")
        her13conv_fn = read_conf(conf,"her13conv_fn")

        soma_5_to_1_3K_fh = project_path*annot_sub_path*her13conv_fn
        
        soma_5_to_1_3K_df = CSV.read(soma_5_to_1_3K_fh,DataFrame)
        soma_5_to_1_3K_dict = Dict()
        for i in 1:size(soma_5_to_1_3K_df)[1]
            soma_5_to_1_3K_dict[soma_5_to_1_3K_df.heritage_name[i]] = soma_5_to_1_3K_df.soma_id[i]
        end

        ## CREATE var_key and check for repeated var names
        var_key = Dict()
        not_unique = []
        for name in prot_names
            var_key[name] = name
        end
        for name in prot_names
            if startswith(name,r"[A-z]_")
                var_key[name] = name[3:end]
                if haskey(soma_5_to_1_3K_dict,name[3:end]) ## if conversion dict has name, then use that value (e.g. without A/B)
                    #println("Convert $(name) to 1.3K")
                    var_key[name] = soma_5_to_1_3K_dict[name[3:end]]
                end
            end
        end
        for name in met_names
            idxs = []
            for i in 1:length(ogm_key_ami_df[:,Symbol(study*study_suffix*"_name")])
                if !ismissing(ogm_key_ami_df[i,Symbol(study*study_suffix*"_name")])
                    if ogm_key_ami_df[i,Symbol(study*study_suffix*"_name")] == name || ogm_key_ami_df[i,Symbol(study*study_suffix*"_name")]*" " == name
                        push!(idxs,i)
                    end
                end
            end
            if length(idxs) == 1
                var_key[name] = ogm_key_ami_df.ogm_id[idxs[1]]
            elseif length(idxs) > 1
                println("VAR $name is not unique in ogm file!!!")
                push!(not_unique,name)
            else
                var_key[name] = ""
            end
        end
        for name in met_names
            idxs = []
            for i in 1:length(ogm_key_hil_df[:,Symbol(study*study_suffix*"_name")])
                if !ismissing(ogm_key_hil_df[i,Symbol(study*study_suffix*"_name")])
                    if ogm_key_hil_df[i,Symbol(study*study_suffix*"_name")] == name || ogm_key_hil_df[i,Symbol(study*study_suffix*"_name")]*" " == name
                        push!(idxs,i)
                    end
                end
            end
            if length(idxs) == 1
                var_key[name] = ogm_key_hil_df.ogm_id[idxs[1]]
            elseif length(idxs) > 1
                println("VAR $name is not unique in ogm file!!!")
                push!(not_unique,name)
            elseif !haskey(var_key,name)
                var_key[name] = ""
            end
        end

        # Compute correlations

        cor_filt_cat_df = DataFrame()
        for study_file in sort(study_files_arr,rev = true)
            study_file_fh = project_path*project_sub_path*study_file
            println(study_file_fh)
            study_join_df = CSV.read(study_file_fh,DataFrame)

            study_file_cor_df = find_correlations(prot_names,met_names,study_join_df,var_key)

            # determine variable annotations
            suffix = study_file
            prefix = study*"_join_age_sex"
            suffix = replace(suffix, prefix => "")
            suffix = replace(suffix, "_norm.csv" => "")

            study_prefix_names = study .* "_" .* names(study_file_cor_df)
            rename!(study_file_cor_df, study_prefix_names)

            if suffix != "_"
                study_suffix_names = names(study_file_cor_df) .* suffix
                rename!(study_file_cor_df, study_suffix_names)
            end            

            cor_fn = split(study_file,".")[1] * "_cor.csv"
            cor_fh = project_path*project_sub_out_path*cor_fn
            #CSV.write(cor_fh,study_file_cor_df)

            cor_filt_cat_df = hcat(cor_filt_cat_df,study_file_cor_df,makeunique = true)
        end

        cor_filt_cat_df = filter(row -> row[Symbol(study*"_var1_alt")] != "" && row[Symbol(study*"_var2_alt")] != "",cor_filt_cat_df)

        cor_filt_cat_df[:,:cor_id] .= ""
        for i in 1:size(cor_filt_cat_df)[1]
            cor_id_arr = sort([cor_filt_cat_df[i,Symbol(study*"_var1_alt")],cor_filt_cat_df[i,Symbol(study*"_var2_alt")]])
            cor_filt_cat_df.cor_id[i] = join(cor_id_arr)
        end

        study_cor_filt_fh = project_path*project_sub_out_path*study*"_cor_filt.csv"
        CSV.write(study_cor_filt_fh,cor_filt_cat_df)

    end

end

studies_arr = ["jhs","mesa","her"]
#studies_arr = ["mesa"]
pm_correlations(studies_arr)