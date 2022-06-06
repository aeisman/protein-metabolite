Distributed.@everywhere using RCall,CSV,DataFrames,GZip,Random
Distributed.@everywhere include("gwas_tools.jl")

function call_gwas_mr(root_path,use_corr,fdr_level,tss_tol,study)
    # Load ogm meta anno key
    anno_path = "$(root_path)/data/annotations/"
    ogm = CSV.read(anno_path*"ogm_meta_key_anno.csv",DataFrame)
    ogm_study_dict = Dict()

    # Load study specific correlations file for LIML MR method
    cor_path = "$(root_path)/data/correlations/"
    if study == "jhs"
        cor_all = CSV.read(cor_path*"jhs_cor_filt.csv",DataFrame);
        rename!(cor_all,:jhs_rho => :rho);
        cor = select(cor_all,[:cor_id,:rho]);

        filter!(row -> !ismissing(row.jhs_id), ogm)
        for i in 1:size(ogm)[1]
            ogm_study_dict[ogm.jhs_id[i]] = ogm.ogm_id[i]
        end
    elseif study == "mesa"
        cor_all = CSV.read(cor_path*"mesa_cor_filt.csv",DataFrame);
        rename!(cor_all,:mesa_rho => :rho);
        cor = select(cor_all,[:cor_id,:rho]);

        filter!(row -> !ismissing(row.mesa_id), ogm)
        for i in 1:size(ogm)[1]
            ogm_study_dict[ogm.mesa_id[i]] = ogm.ogm_id[i]
        end
    elseif study == "heritage"
        cor_all = CSV.read(cor_path*"her_cor_filt.csv",DataFrame);
        rename!(cor_all,:her_rho => :rho);
        cor = select(cor_all,[:cor_id,:rho]);

        filter!(row -> !ismissing(row.heritage_id), ogm)
        for i in 1:size(ogm)[1]
            ogm_study_dict[ogm.heritage_id[i]] = ogm.ogm_id[i]
        end
    end

    #==
    if use_corr == 1
        ld_dict = gwas_ld_dict(root_path,fdr_level,tss_tol*2,study)
    else
        ld_dict = 0
    end
    ==#

    #using RCall,CSV,DataFrames,GZip
    #include("omics_tools.jl")

    #gwasjoin_df = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/data/gwasjoin/SL004768_SL005250_0.1.shared.gwasjoin") |> DataFrame;
    #exp_name = "SL004768"
    #out_name = "SL005250"

    files_arr = readdir("$(root_path)/data/gwas_protein_out", sort = false);
    #shared_done = readdir("$(root_path)/data/gwas_protein_out/gwasjoin", sort = false);
    gwasjoin = files_arr[endswith.(files_arr,".gwasjoin")];

    ldr_files_arr = readdir("$(root_path)/data/gwas_protein_out/ldr", sort = false);
    ldr_matrix_files = ldr_files_arr[endswith.(ldr_files_arr,".ldr_snps_plink.ld")];

    gwasjoin_keep = []
    for file in gwasjoin
        gwasjoin_fh = "$(root_path)/data/gwas_protein_out/$(file)"
        file_arr = split(file,"_")
        if use_corr == 1
            try
                ldr_matrix_file = ldr_matrix_files[contains.(ldr_matrix_files,file_arr[1])][1]
                if !isfile(gwasjoin_fh*".mr") && isfile("$(root_path)/data/gwas_protein_out/ldr/$(ldr_matrix_file)")
                    push!(gwasjoin_keep,file)
                end
            catch
            end
        else
            push!(gwasjoin_keep,file)
        end
    end
    gwasjoin = gwasjoin_keep

    n_to_do = length(gwasjoin)
    println("$(n_to_do) to do!")

    #shuffle!(gwasjoin)
    
    @Distributed.distributed vcat for file in gwasjoin
        println(file)
        gwasjoin_fh = "$(root_path)/data/gwas_protein_out/$(file)"
        if !isfile(gwasjoin_fh*".mr")
            file_arr = split(file,"_")
            file_arr[2] = split(file_arr[2],".")[1]
            if file_arr[1] != file_arr[2]
                gwasjoin_df = CSV.read(gwasjoin_fh,DataFrame);

                ld_r = []
                if use_corr == 1
                    try
                        ldr_matrix_file = ldr_matrix_files[contains.(ldr_matrix_files,file_arr[1])][1]
                        ld_r = CSV.read("$(root_path)/data/gwas_protein_out/ldr/$(ldr_matrix_file)",DataFrame,header=false)
                    catch
                        @warn "No LD file for $(file_arr)[1]"
                    end
                end
                #if size(gwasjoin_df)[1] >= 3
                if ld_r != [] || use_corr != 1
                    println("Trying $(gwasjoin_fh)...")
                    try
                        cor_id_f1 = file_arr[1]
                        cor_id_f2 = file_arr[2]
                        if haskey(ogm_study_dict, cor_id_f1)
                            cor_id_f1 = ogm_study_dict[cor_id_f1]
                        end
                        if haskey(ogm_study_dict, cor_id_f2)
                            cor_id_f2 = ogm_study_dict[cor_id_f2]
                        end

                        cor_id = join(sort([cor_id_f1,cor_id_f2]))
                        psi = 999
                        try
                            psi = filter(row -> row.cor_id == cor_id, cor).rho[1]
                        catch
                            @warn("No Correlation Available for LIML Method!")
                        end

                        #println("Call gwas_mr")
                        #mrv_row = gwas_mr(gwasjoin_df,file_arr[1],file_arr[2],use_corr,ld_r)
                        println("Starting MR for $cor_id with PSI $psi...")
                        mrv_row = gwas_mr(gwasjoin_df,file_arr[1],file_arr[2],use_corr,ld_r,psi) # add psi for LIML 4/25/22
                        println("Finished MR...")
                        #println("Finished gwas_mr")
                        #mrv_row = gwas_mr(gwasjoin_df,file_arr[1],file_arr[2],use_corr,ld_dict)
                        #println("Write mr result")
                        CSV.write(gwasjoin_fh*".mr",mrv_row, writeheader=false)
                        #println("MR result written!")
                    catch
                        @warn "$(gwasjoin_fh) failed!"
                    end
                #end
                end
            end
        end
    end

end

#root_path = "/home/aaron_eisman"
#root_path = "/home/aaroneisman"
#call_gwas_mr(root_path)