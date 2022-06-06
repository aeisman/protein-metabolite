using CSV, DataFrames

data_path = "/home/aaroneisman/data/"
results_sub_folder = "20220425/"
path = data_path*results_sub_folder

gsutil_str = "gsutil cp -r gs://jhs_data_topmed/eisman/20220425 $(data_path)"
run(`bash -c "$(gsutil_str)"`)

# get names of .tar.gz files in that folder
files = readdir(path)
to_extract = files[endswith.(files,".tar.gz")]
for file in to_extract
    println(file)
    pigz_str = "unpigz < $(data_path)$(results_sub_folder)$(file) | tar -xvC ../data/20220425"
    run(`bash -c "$(pigz_str)"`)
end

her_key_df = CSV.read("../data/heritage_13soma_conv.csv",DataFrame)
her_key_dict = Dict()
for i in 1:size(her_key_df)[1]
    her_key_dict[her_key_df.heritage_name[i]] = her_key_df.soma_id[i]
end

function get_instruments(path, folder, her_key_dict)

    instr_keep_dict = Dict()
    instr_select_dict = Dict()

    files_all = readdir(path*folder)
    files_keep = files_all[contains.(files_all,"keep")]

    for file in files_keep
        if startswith(file,"B_")
            study_id = split(file,"_")[2]
        else
            study_id = split(file,"_")[1]
        end
        t_study_id_dict = Dict()
        #println(study_id)

        fh = path*folder*"/"*file

        t_df = CSV.read(fh,DataFrame)
        #for i in 1:size(t_df)[1]
        #    t_study_id_dict[t_df.SNP[i]] = t_df[i,:]
        #end
        t_study_id_dict = t_df
        instr_keep_dict[study_id] = t_study_id_dict
    end

    files_select = readdir(path*folder*"/ldr/")
    files_select = files_select[endswith.(files_select,".ldr")]

    for file in files_select
        #prot = match(r"SL[0-9]*",file).match
        prot = match(r"SL[0-9]*[AB]?(?=_)",file).match;
        fh = path*folder*"/ldr/"*file;
        t_df = CSV.read(fh,DataFrame);
        instr_select_dict[prot] = Set(t_df.SNP);
        #for line in readlines(fh)
        #    if haskey(instr_select_dict,prot)
        #        push!(instr_select_dict[prot],line)
        #    else
        #        instr_select_dict[prot] = Set([line])
        #    end
        #end
    end

    return (instr_select_dict, instr_keep_dict)

end

folders = []
for file in to_extract
    t_folder = file[1:end-7]
    push!(folders,t_folder)
end

for folder in folders
    println()
    println("### STARTING $(folder) ###")
    print("Generating instrument dictionaries...")
    (t_instr_select_dict, t_instr_keep_dict) = get_instruments(path, folder, her_key_dict)
    println("Finished!!!")
    println()

    t_prots = collect(keys(t_instr_select_dict))
    t_all = collect(keys(t_instr_keep_dict))

    println("Generating MR instrument table...")
    mr_instrument_df = DataFrame()
    for prot in t_prots
        println(prot)
        t_prot_snps = t_instr_select_dict[prot]
        if length(t_prot_snps)[1] > 0
            t_prot_select_gwas = filter(row -> in(row.SNP, t_prot_snps), t_instr_keep_dict[prot])
            for outcome in t_all
                if prot != outcome
                    try
                        t_outcome_select_gwas = filter(row -> in(row.SNP, t_prot_snps), t_instr_keep_dict[outcome])
                        t_mr_instrument_df = leftjoin(t_prot_select_gwas,t_outcome_select_gwas,on = :SNP, makeunique = true)
                        t_mr_instrument_df[!,:exposure] .= prot
                        t_mr_instrument_df[!,:outcome] .= outcome
                        #mr_instrument_df = vcat(mr_instrument_df,t_mr_instrument_df)
                        append!(mr_instrument_df,t_mr_instrument_df)
                    catch
                        @warn("Outcome $(outcome) failed!!!")
                    end
                end
            end
        end
    end

    println("Writing MR instrument table to file")
    CSV.write(path*folder*"_mr_instrument_df.csv",mr_instrument_df)
end
