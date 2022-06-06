### Determine which heritage proteins with two soma_ids (A and B) line up with the 1.3K platform

using DataFrames,CSV

jhs = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/data/annotations/heritage_1.3K_conversion/JHS_protein_ID_key.csv",DataFrame)
jhs.SeqId = replace.(jhs.SeqId,r"_.*" => "") ### remove _# part of jhs SeqId
#RegexMatch.(r".*(?=_)",jhs.SeqId)

her = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/data/annotations/heritage_1.3K_conversion/HERITAGE Somalogic analytes.csv",DataFrame)


jhs_seq_set = Set(jhs.SeqId)
her_seq_set = Set(her.SeqId)
intersect(jhs_seq_set,her_seq_set) # 1159

jhs_soma_ids = Set(jhs.SomaId)
her_all_soma_ids = Set(her.SomaId)
intersect(her_all_soma_ids,jhs_soma_ids) # 1297


### filter heritage by jhs soma soma_ids
her_filt = filter(row -> in(row.SomaId, jhs_soma_ids),her)
her_filt2 = filter(row -> length(row.SASNAME) == 9, her_filt) ## names with an A or a B at the end


her_filt2.SeqId

her_filt2[!,:SeqId_in_jhs] .= 9
for i in 1:size(her_filt2)[1]
    her_filt2.SeqId_in_jhs[i] = in(her_filt2.SeqId[i],jhs_seq_set)
end

her_somaid_w_jhs_match = Set(filter(row -> row.SeqId_in_jhs == 1, her_filt2).SomaId)

filter(row -> !in(row.SomaId,her_somaid_w_jhs_match),her_filt2)

# manual
manual_selects = Set(["SL003722A","SL004580A","SL004864B","SL005222B","SL007429A"])
her_ab_selects = filter(row -> row.SeqId_in_jhs == 1 || in(row.SASNAME,manual_selects),her_filt2)
heritage_key = DataFrame(soma_id = her_ab_selects.SomaId, heritage_name = her_ab_selects.SASNAME)
CSV.write("heritage_13soma_conv.csv",heritage_key)




#== There are some cases where there has been a one for one seqid swap in the move from 1.3K to 5k
    study = "heritage"
    data_folder = "./"
    on_google = 1
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

somaid_arr_keep = ["SL003722A","SL004580A","SL004580B","SL004864A","SL004864B","SL005222A","SL005222B","SL007429A","SL007429B"]
somaid_arr_keep = ["SL007429A","SL007429B"]
for i in 1:length(somaid_arr_keep)
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
    catch
    end
end



function min_p(gwas_fh)
    gwas_io = GZip.open(gwas_fh)
    p_min = 1.0
    for line in eachline(gwas_io)
        try
            p = parse(Float64,split(line,"\t")[10])
            
            if p < p_min
                p_min = p

            end
        catch
            println("skipped line")
        end
    end
    close(gwas_io)
    return(p_min)
end

files = readdir(".")
for file in files
    #if startswith(file,"B_")
    if contains(file,"7429")
        println(file)
        p = min_p(file)
        println(p)
        println("\n\n")
    end
end


B_SL003722A.txt.fastgwa.gz
skipped line
2.34748e-16


B_SL004580A.txt.fastgwa.gz
skipped line
6.5807e-9


B_SL004580B.txt.fastgwa.gz
skipped line
1.12342e-8


B_SL004864A.txt.fastgwa.gz
skipped line
4.83183e-7


B_SL004864B.txt.fastgwa.gz
skipped line
5.86702e-11


B_SL005222A.txt.fastgwa.gz
skipped line
1.33656e-8


B_SL005222B.txt.fastgwa.gz
skipped line
4.64038e-9


B_SL007429A.txt.fastgwa.gz
skipped line
1.79284e-10


B_SL007429B.txt.fastgwa.gz
skipped line
3.85901e-7

=#