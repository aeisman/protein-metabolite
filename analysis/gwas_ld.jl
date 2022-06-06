Distributed.@everywhere using RCall,CSV,DataFrames,GZip,Random
Distributed.@everywhere include("gwas_tools.jl")

function call_gwas_ld_reduce(root_path,ld_thresh,fdr_level,gwas_file_post,tss_tol,maf,study)
    if study == "jhs"
        ## JHS
        bfile_fh = "/home/aaroneisman/data/genotypes/jhs_proteins"
    elseif study == "mesa"
        ## MESA
        bfile_fh = "/home/aaroneisman/data/genotypes/mesa_metab"
    elseif study == "heritage"
        bfile_fh = "/home/aaroneisman/data/genotypes/heritage_merge"
    else
        @error("Study $(study) not recognized!")
    end

    plink_ld_window = tss_tol / 1000
    #==
    #load SNP list
    #root_path = "/Users/aaroneisman/projects/benson"
    #root_path = "/home/aaron_eisman"

    #JHS
    all_fdr_df = CSV.read("$(root_path)/data/gwas_protein_out/jhs_prot_adjust_age_sex_pcs_prot_all_fdr_$(fdr_level).txt",DataFrame)

    #MESA
    #all_fdr_df = CSV.read("$(root_path)/data/gwas_protein_out/all_fdr_$(fdr_level).txt",DataFrame);
    ##all_shared_df = CSV.read("$(root_path)/data/google_cloud_out/gwas_protein_out_20200913/all.shared.csv") |> DataFrame

    all_snp_set = Set(all_fdr_df.SNP);

    #ld_fh = "/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/data/gwas_protein/jhs_geno_plink.ld"
    #ld_fh = "/home/aaron_eisman/data/plink.ld"
    ld_fh = "$(root_path)/data/plink.001.ld"
    ld_io = open(ld_fh,"r")
    ld_dict = Dict()
    for line in readlines(ld_io)
        #println(line)
        line_arr = split(line,r"\s+")
        temp_snp_arr = [line_arr[4],line_arr[7]]
        temp_snp_set = Set(temp_snp_arr)
        #sort!(snp_arr)
        #snp_str = snp_arr[1]*snp_arr[2]
        #println(snp_str)

        if in(temp_snp_arr[1],all_snp_set) && in(temp_snp_arr[2],all_snp_set)
            #println(line)
            ld_dict[temp_snp_set] = parse(Float64,line_arr[8])
        end
    end
    close(ld_io)
    ==#

    #=
    1) get list of .shared files
    2) make ldr folder
    3) ld_reduce each file (distributed) and move to ldr folder
    =#
    files_arr = readdir("$(root_path)/data/gwas_protein_out", sort = false);
    #shared_done = readdir("$(root_path)/data/gwas_protein_out/ldr", sort = false);
    fdr_files = files_arr[endswith.(files_arr,"_fdr_$(fdr_level)_maf_$(maf)$(gwas_file_post)")];

    gene_locs_file = "$(root_path)/data/gene_locs.csv"

    gene_locs_df = CSV.read(gene_locs_file,DataFrame)
    gene_locs_dict = Dict()

    for row in eachrow(gene_locs_df)
        gene_locs_dict[row.somaid] = [row.chromosome,row.tss]        
    end

    try
        mkdir("$(root_path)/data/gwas_protein_out/ldr")
    catch
        @warn "Folder already exists!"
    end

    println("Checking for overlapping cis regions and already completed files")
    fdr_keep = []
    for file in fdr_files
        file_arr = split(file,"_")

        if (match(r"^M",file_arr[1])) != nothing && (match(r"^M",file_arr[2])) != nothing
            sl1_chr = gene_locs_dict[file_arr[1]][1]
            sl2_chr = gene_locs_dict[file_arr[2]][1]
            sl1_tss = gene_locs_dict[file_arr[1]][2]
            sl2_tss = gene_locs_dict[file_arr[2]][2]

            if (sl1_chr == sl2_chr) && (abs(sl1_tss-sl2_tss) < tss_tol)
                #println("$(file_arr[1]) and $(file_arr[2]) have the same TSS, skipping!")
                #mv("$(root_path)/data/gwas_protein_out/$(file)","$(root_path)/data/gwas_protein_out/$(file).skipped")
                
            else
                if !isfile("$(root_path)/data/gwas_protein_out/ldr/$(file).ldr")
                    push!(fdr_keep,file)
                end
            end
        else
            if !isfile("$(root_path)/data/gwas_protein_out/ldr/$(file).ldr")
                push!(fdr_keep,file)
            end
        end
    end

    shuffle!(fdr_keep)
    

    println("Starting LD Reduce")
    if length(fdr_keep) > 0
        @Distributed.distributed vcat for file in fdr_keep ## consider breaking up into groups of 100 to reduce total workers generated
        #@Distributed.distributed vcat for i in 1:100:length(shared) ## consider breaking up into groups of 100 to reduce total workers generated
        #    shared_chunk = shared[i:min(i+99,length(shared))]
        #    for file in shared_chunk
            #for i in 1:1

                gwas_fh = "$(root_path)/data/gwas_protein_out/$(file)"
                #println("Reducing $(gwas_fh)")
                if !isfile("$(root_path)/data/gwas_protein_out/ldr/$(file).ldr")
                    println("Start LD for: $(gwas_fh)")
                    try
                        keep_df = CSV.read(gwas_fh,DataFrame)

                        if size(keep_df)[1] > 1
                            ###snp_col = 2
                            ###p_col = 10
                            ###gwas_ld_df = gwas_ld_clump(gwas_fh,snp_col,p_col,ld_dict,ld_thresh)

                            #use sed to convert to tab delimited file for plink
                            #sed -i 's/,/\t/g' /home/aaroneisman/data/gwas_protein_out_fdr_0.05_ld_0.01_tss_1000000/jhs_prot_adjust_age_sex_pcs_prot_SL005574_fdr_0.05_maf_0.01_cis.txt >/home/aaroneisman/data/gwas_protein_out_fdr_0.05_ld_0.01_tss_1000000/jhs_prot_adjust_age_sex_pcs_prot_SL005574_fdr_0.05_maf_0.01_cis.txt.tab
                            sed1_str = "sed -i 's/,/\t/g' $(gwas_fh)"
                            run(`bash -c "$(sed1_str)"`)

                            #run plink
                            #plink --bfile /home/aaroneisman/data/genotypes/jhs_proteins --clump /home/aaroneisman/data/gwas_protein_out_fdr_0.05_ld_0.01_tss_1000000/jhs_prot_adjust_age_sex_pcs_prot_SL005574_fdr_0.05_maf_0.01_cis.txt --clump-r2 0.001 --clump-kb 1000 --memory 8000 --out /home/aaroneisman/data/gwas_protein_out_fdr_0.05_ld_0.01_tss_1000000/jhs_prot_adjust_age_sex_pcs_prot_SL005574_fdr_0.05_maf_0.01_cis.txt
                            plink_str = "plink --bfile $(bfile_fh) --clump $(gwas_fh) --clump-r2 $(ld_thresh) --clump-kb $(plink_ld_window) --memory 8000 --keep-allele-order --out $(gwas_fh)"
                            run(`bash -c "$(plink_str)"`)

                            #use sed to convert plink output to csv
                            #sed -i -E 's/(\s\s*+)/,/g' /home/aaroneisman/data/gwas_protein_out_fdr_0.05_ld_0.01_tss_1000000/jhs_prot_adjust_age_sex_pcs_prot_SL005574_fdr_0.05_maf_0.01_cis.txt.clumped
                            sed2_str = "sed -i -E 's/(,)/|/g' $(gwas_fh).clumped"
                            run(`bash -c "$(sed2_str)"`)
                            sed3_str = "sed -i -E 's/(\\s\\s*+)/,/g' $(gwas_fh).clumped"
                            run(`bash -c "$(sed3_str)"`)

                            #read plink output and reduce gwas file
                            clump_df = CSV.read("$(gwas_fh).clumped",DataFrame)
                            clump_snps = Set(clump_df.SNP)
                            filter!(row -> row.SNP in clump_snps, keep_df)
                        end

                        # write LD reduced file
                        CSV.write("$(root_path)/data/gwas_protein_out/ldr/$(file).ldr",keep_df)
                    catch
                        println("LD Clump of $(gwas_fh) failed!!!")
                        @warn "LD Clump of $(gwas_fh) failed!!!"
                    end
                    println("Finished LD for: $(gwas_fh)")
                    #mv("$(gwas_fh).ldr","$(root_path)/data/gwas_protein_out/ldr/$(file).ldr")
                end
            #end
        end
    end

    cat_str = "cat"

    files_arr = readdir("$(root_path)/data/gwas_protein_out/ldr", sort = false);
    ldr_files = files_arr[endswith.(files_arr,".ldr")];

    for file in ldr_files
        gwas_keep_fh = "$(root_path)/data/gwas_protein_out/ldr/"*file
        cat_str *= " $(gwas_keep_fh)"
    end

    cat_str *= " >"*"$(root_path)/data/gwas_protein_out/ldr/"*"all_ldr_$(fdr_level).txt"

    ### execute cat command
    run(`bash -c "$(cat_str)"`)

    ## Generate Correlation Matrices for all LD files
    files_arr = readdir("$(root_path)/data/gwas_protein_out/ldr", sort = false)
    ldr_files = files_arr[endswith.(files_arr,".ldr")]
    @Distributed.distributed vcat for file in ldr_files
        try
            if !isfile("$(root_path)/data/gwas_protein_out/ldr/$(file)_snps_plink.ld")
                gwas_ld_r(root_path,file,study)
            end
        catch
            @warn("Failed gwas_ld_r for $(file)")
        end
    end

end

#root_path = "/home/aaron_eisman"
#root_path = "/home/aaroneisman"
#ld_thresh = 0.01
#fdr_level = 0.05
#gwas_file_post = "_cis.txt"
#tss_tol = 1000000
#call_gwas_ld_reduce(root_path,ld_thresh,fdr_level,gwas_file_post,tss_tol)