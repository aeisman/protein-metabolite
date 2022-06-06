Distributed.@everywhere using RCall,CSV,DataFrames,Random
Distributed.@everywhere include("gwas_tools.jl")

function call_gwas_maf_filter(somaid_arr,maf,root_path,study)
    data_folder = "$(root_path)/data/gwas_protein_out/"

    if study == "jhs"
        #JHS
        gwas_file_pre = "jhs_prot_adjust_age_sex_pcs_prot_"
    elseif study == "mesa"
        #MESA
        gwas_file_pre = ""
    elseif study == "heritage"
        gwas_file_pre = "B_"
    else
        @error("Study $(study) not recognized!")
    end

    gwas_file_post = "_cis.txt"

    gwas_af1_col = 7

    Distributed.@distributed vcat for somaid in somaid_arr
        println("Starting MAF Filter for $(somaid)")
        gwas_fh = data_folder*gwas_file_pre*somaid*gwas_file_post
        if stat(gwas_fh).size > 0
            gwas_maf_fh = data_folder*gwas_file_pre*somaid*"_maf_$(maf)"*gwas_file_post
            gwas_maf_filter(gwas_fh,gwas_maf_fh,gwas_af1_col,maf)
            println("Finished MAF Filter for $(somaid)")
        else
            println("Failed MAF Filter for $(somaid), GWAS file was size 0!!!")
        end
    end

end