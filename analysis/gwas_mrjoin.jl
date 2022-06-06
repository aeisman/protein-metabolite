Distributed.@everywhere using RCall,CSV,DataFrames,Random
Distributed.@everywhere include("gwas_tools.jl")


function call_gwas_innerjoin_mr(root_path,fdr_level,maf,study)

    gwas_file_post = "_cis.txt"

    files_arr = readdir("$(root_path)/data/gwas_protein_out", sort = false);
    exp_files = files_arr[endswith.(files_arr,"_fdr_$(fdr_level)_maf_$(maf)$(gwas_file_post)")];
    out_files = files_arr[endswith.(files_arr,"_keep.fastgwa.txt")];

    shuffle!(exp_files)
    shuffle!(out_files)

    n_to_do = length(exp_files) * length(out_files)
    println("$n_to_do to do!!!")

    Distributed.@distributed vcat for exp_file in exp_files
        exp_file_arr = split(exp_file,"_")
        
        if study == "jhs"
            ##JHS
            exp_file_short = exp_file_arr[8]
        elseif study == "mesa"
            ##MESA
            exp_file_short = exp_file_arr[1] ## location of somiad in filename
        elseif study == "heritage"
            ##HERITAGE
            exp_file_short = exp_file_arr[2]
        else
            @error("Study $(study) not recognized!")
        end

        gwas_fh_a = "$(root_path)/data/gwas_protein_out/$(exp_file)"

        println("Starting $(exp_file_short)...")
        for out_file in out_files
            println("Out file:"*out_file)
            out_file_arr = split(out_file,"_")

            if study == "jhs"
                ##JHS
                if length(out_file_arr) > 3
                    #JHS proteins
                    out_file_short = out_file_arr[8]
                else
                    #JHS Metabolites
                    out_file_short = out_file_arr[1]
                end
            elseif study == "mesa"
                ##MESA
                out_file_short = out_file_arr[1]
            elseif study == "heritage"
                ##HERITAGE
                ### problem for protein-protein
                if out_file_arr[1] == "B"
                    out_file_short = out_file_arr[2]
                else    
                    out_file_short = out_file_arr[1]
                end
            else
                @error("Study $(study) not recognized!")
            end

            #println(out_file_short)

            gwas_fh_b = "$(root_path)/data/gwas_protein_out/$(out_file)"
            out_fh = "$(root_path)/data/gwas_protein_out/$(exp_file_short)_$(out_file_short).gwasjoin"
            snp_col = 2
            try
                gwas_innerjoin(gwas_fh_a,snp_col,gwas_fh_b,snp_col,out_fh)
                gwasjoin_df = CSV.read(out_fh,DataFrame)
                
                ldr_df = CSV.read("$(root_path)/data/gwas_protein_out/ldr/$(exp_file).ldr",DataFrame)
                ldr_snps = Set(ldr_df.SNP)
                filter!(row -> row.SNP in ldr_snps, gwasjoin_df)

                CSV.write(out_fh, gwasjoin_df)

            catch
                @warn "$(exp_file_short) join with $(out_file_short) failed!!!"
            end

        end
        println("Finished $(exp_file_short)!!!")
    end

end

#root_path = "/home/aaroneisman"
#root_path = "/home/aaron_eisman"
#fdr_level = 0.05

#call_gwas_innerjoin_mr(root_path,fdr_level)