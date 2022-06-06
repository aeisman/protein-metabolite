Distributed.@everywhere using RCall,CSV,DataFrames,Random
Distributed.@everywhere include("gwas_tools.jl")

### test arguments
function call_gwas_fdr(somaid_arr,fdr_level,root_path,maf,study)
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

    gwas_file_post = "_maf_$(maf)_cis.txt"


    #fdr_level = 0.1

    #somaid_arr = ["SL005699","SL012707","SL000020"]
    #somaid_arr = ["SL000020"]

    cat_str = "cat"

    shuffle!(somaid_arr)

    Distributed.@distributed vcat for somaid in somaid_arr
        println("Starting FDR for $(somaid)")
        gwas_fh = data_folder*gwas_file_pre*somaid*gwas_file_post
        if stat(gwas_fh).size > 0
            gwas_keep_fh = data_folder*gwas_file_pre*somaid*"_fdr_$(fdr_level)"*gwas_file_post
            if fdr_level == "bonf"
                ### if FDR level is "bonf" then do a bonferoni significance level instead of fdr method
                sig = 0.05 ## will be divided by number of p values being tested
                gwas_bonf(gwas_fh, gwas_keep_fh, sig)
            else
                gwas_fdr(gwas_fh,gwas_keep_fh,fdr_level)
            end
            println("Finished FDR for $(somaid)")
        else
            println("Failed FDR for $(somaid), GWAS file was size 0!!!")
        end
    end

    files_arr = readdir("$(root_path)/data/gwas_protein_out", sort = false);
    fdr_files = files_arr[endswith.(files_arr,"_fdr_$(fdr_level)$(gwas_file_post)")];

    for file in fdr_files
        gwas_keep_fh = data_folder*file
        cat_str *= " $(gwas_keep_fh)"
    end

    cat_str *= " >"*data_folder*gwas_file_pre*"all_fdr_$(fdr_level).txt"

    ### execute cat command
    run(`bash -c "$(cat_str)"`)
    #@rput(cat_str)
    #R"""
    #system(cat_str)
    #"""
end

#call_gwas_fdr(somaid_arr,fdr_level,root_path)