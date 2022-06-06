Distributed.@everywhere using RCall,CSV,DataFrames,Random
Distributed.@everywhere include("gwas_tools.jl")

function call_gwas_cis_regions(somaid_arr,root_path,on_google,tss_tol,study)
    ### test arguments
    # SL000020 Apolipoprotein B
    #p = "SL000324"
    #p = "SL005699" # Apolipoprotein L1 (variant in L1 that ... classic positve selection)
    #B or E
    #p = "SL012707" # PCSK9

    #somaid_arr = ["SL005699","SL012707","SL000020","SL000276"]
    #somaid_arr = ["SL000020"]


    #n_cis_regions = length(cis_regions_set)
    #gwas_file = "../data/dan/jhs_prot_adjust_age_sex_pcs_prot_SL005699.fastgwa.gz"

    # Load gene locs file
    gene_locs_file = "$(root_path)/data/gene_locs.csv"
    gene_locs_df = CSV.read(gene_locs_file, DataFrame)

    # Create set of somaid with known gene locs
    gene_locs_set = Set(gene_locs_df.somaid)
    
    # Create set of somaid's from somaid input array
    somaid_arr_set = Set(somaid_arr)
    # Determine which of the somaid input array have known gene locs
    #somaid_arr = collect(intersect(somaid_arr_set,gene_locs_set)) ## creates problems with A/B proteins

    prot_have_set = Set()
    for somaid in somaid_arr
        her_prot = somaid ### changed to include AB proteins on 1/12/2022
        prot = match(r"SL[0-9]*",somaid).match
        #push!(prot_have_set,chop(prot))
        if in(prot,gene_locs_set)
            push!(prot_have_set,her_prot)
        end
    end

    somaid_arr = collect(prot_have_set)
    
    # Define variables about data location and data files
    data_folder = "$(root_path)/data/gwas_protein/"
    out_folder = "$(root_path)/data/gwas_protein_out/"
    
    if study == "jhs"
        #JHS
        gwas_file_pre = "jhs_prot_adjust_age_sex_pcs_prot_"
        gwas_file_post = ".fastgwa.gz"
    elseif study == "mesa"
        #MESA
        gwas_file_pre = ""
        gwas_file_post = "_1_resid.txt_baseline.fastgwa.gz"
    elseif study == "heritage"
        #HERITAGE
        gwas_file_pre = "B_"
        gwas_file_post = ".txt.fastgwa.gz"
    else
        @error("Study $(study) not recognized!")
    end

    gwas_chr_col = 1
    gwas_pos_col = 3
    #tss_tol = 1000000 # transcription start site tolerance for "cis" region

    somaid_arr_keep = []
    for somaid in somaid_arr
        gwas_cis_fh = out_folder*gwas_file_pre*somaid*"_cis.txt"
        if !isfile(gwas_cis_fh)
            push!(somaid_arr_keep,somaid)
        end
    end

    shuffle!(somaid_arr_keep)

    if length(somaid_arr_keep) > 0
        Distributed.@distributed vcat for i in 1:length(somaid_arr_keep)
            gwas_file = "$(gwas_file_pre)$(somaid_arr_keep[i])$(gwas_file_post)"

            try
                println("Starting cis-region extraction for $(somaid_arr_keep[i])")
                if on_google == 1
                    # if running on google, copy a temporary data file to work with from bucket
                    if study == "jhs"
                        #JHS
                        run(`gsutil cp gs://jhs_data_topmed/jhs_gwas/jhs_SOMA_gwas/fastgwa_with_pcs/raw_results/$(somaid_arr_keep[i])$(gwas_file_post) $(data_folder)$(gwas_file)`)
                    elseif study == "mesa"
                        #MESA
                        run(`gsutil cp gs://mesa-bucket/mesa_gwas/mesa_SOMA_gwas/fastgwa/raw_results/$(somaid_arr_keep[i])$(gwas_file_post) $(data_folder)$(gwas_file)`)
                    elseif study == "heritage"
                        #HERITAGE
                        run(`gsutil cp gs://heritage-bucket/heritage_gwas/heritage_SOMA_gwas/fastgwa/raw_results/$(gwas_file) $(data_folder)$(gwas_file)`)
                    else
                        @error("Study $(study) not recognized!")
                    end
                end

                somaid = somaid_arr_keep[i]
                gwas_fh = data_folder*gwas_file_pre*somaid*gwas_file_post
                gwas_cis_fh = out_folder*gwas_file_pre*somaid*"_cis.txt"

                tss_col = 3
                gene_end_col = 6
                cis_regions_set = load_cis_regions(gene_locs_file, somaid, 1, 4, tss_col, tss_tol, gene_end_col) # to gene_end
                #cis_regions_set = load_cis_regions(gene_locs_file, somaid, 1, 4, tss_col, tss_tol, tss_col) # don't include gene
                println("Call gwas_cis_regions()...")
                gwas_cis_regions(gwas_fh,gwas_cis_fh,gwas_chr_col,gwas_pos_col,cis_regions_set)

                if on_google == 1
                    # if running on google, delete the temporary data file copied from bucket
                    println("Deleting temporary file...")
                    gwas_file = "$(gwas_file_pre)$(somaid_arr_keep[i])$(gwas_file_post)"
                    run(`rm $(data_folder)$(gwas_file)`)
                end
                println("Finished cis-region extraction for $(somaid_arr_keep[i])!")
            catch
                @warn "Could not get cis regions for $gwas_file"
            end
        end
    end

    return somaid_arr
end