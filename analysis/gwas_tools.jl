using RCall,CSV,DataFrames,GZip,Random

######################
### GWAS FUNCTIONS ###
######################

function load_cis_regions(gene_locs_fh, somaid, id_col, chr_col, tss_col, tss_tol, gene_end_col)
    # Identify all tss associated with a somaid and return the cis regions as a set of arrays

    somaid = somaid[1:8] # chop off A or B designation for gene locs determination

    gene_locs_set = Set() ## array of arrays [chr,lower,upper]
    gene_locs_io = open(gene_locs_fh,"r")
    skip(gene_locs_io,1)
    for line in readlines(gene_locs_io)
        line_arr = split(line,",")
        id = line_arr[id_col]
        if somaid == id
            println(id)
            chr = parse(Int64,line_arr[chr_col])
            tss = parse(Int64,line_arr[tss_col])
            gene_end = parse(Int64,line_arr[gene_end_col])
            end_window = max(gene_end,tss+tss_tol)
            loc_arr = [chr,tss-tss_tol,end_window]
            push!(gene_locs_set,loc_arr)
        end 
    end
    close(gene_locs_io)
    return gene_locs_set
end

function gwas_cis_regions(gwas_fh,gwas_cis_fh,gwas_chr_col,gwas_pos_col,cis_regions_set)
    # generate cis region file from a GWAS file given a set of arrays of cis regions

    ### generate awk command
    awk_str = "gzip -cd '$(gwas_fh)' | awk '"
    n_cis_regions = length(cis_regions_set)
    i = 1
    for region in cis_regions_set
        awk_str *= "(\$$(gwas_chr_col) == $(region[1]) && \$$(gwas_pos_col) > $(region[2]) && \$$(gwas_pos_col) < $(region[3]))"
        
        if i < n_cis_regions
            awk_str *= " || "
        else
            awk_str *= "' >$(gwas_cis_fh)"
        end

        i += 1
    end

    println(awk_str)

    ### execute awk command
    run(`bash -c "$(awk_str)"`)
    #@rput(awk_str)
    #R"""
    #system(awk_str)
    #"""

    ### remove lines with 'nan' // possible source of memory leak
    sed_str = "sed '/nan/d' $(gwas_cis_fh) >$(gwas_cis_fh)_temp"
    mv_str = "mv $(gwas_cis_fh)_temp $(gwas_cis_fh)"
    println("Run sed_str...")
    run(`bash -c "$(sed_str)"`)
    println("Run mv str...")
    run(`bash -c "$(mv_str)"`)
    
    #@rput(sed_str,mv_str)
    #R"""
    #system(sed_str)
    #system(mv_str)
    #"""
end

function gwas_nom_filter(gwas_fh,gwas_nom_fh,nom)
    ### generate awk command
    awk_str = "gzip -cd '$(gwas_fh)' | awk '"
end

function gwas_maf_filter(gwas_fh,gwas_maf_fh,gwas_af1_col,maf)
    awk_str = "awk '"

    if maf < 0 || maf > 0.5
        error("Minimum Allele Frequency (MAF) is not between 0 and 0.5!!!")
    end

    maf_arr = [maf,1-maf]

    awk_str *= "(\$$(gwas_af1_col) > $(maf_arr[1]) && \$$(gwas_af1_col) < $(maf_arr[2]))"
    
    awk_str *= "' $(gwas_fh) >$(gwas_maf_fh)"

    println(awk_str)

    ### execute awk command
    run(`bash -c "$(awk_str)"`)
end

function gwas_fdr(gwas_fh,gwas_keep_fh,fdr_level)
    # Perform FDR analysis on a GWAS file
    
    ### read in gwas file (product of gwas_cis_regions)
    gwas_df = CSV.read(gwas_fh, header=false, DataFrame)
    rename!(gwas_df,["CHR","SNP","POS","A1","A2","N","AF1","BETA","SE","P"]) # assume this format of gwas result file
    gwas_p_df = select(gwas_df, :P) # select p-values to be passed to R for FDR analysis

    gwas_keep = 0
    if size(gwas_p_df)[1] == 1
        println("Only one p-value!")
        #println(gwas_p_df[1])
        if gwas_p_df[1,1] < fdr_level
            gwas_keep = BitArray([1])
        end
    else
        try
            @rput(gwas_p_df,fdr_level,gwas_fh)
            R"""
            library(qvalue)
            graph_fh = paste(gwas_fh,"_pval_graph.png",sep = "")
            csv_fh = paste(gwas_fh,"_pval.csv",sep = "")
            print(graph_fh)
            write.csv(gwas_p_df$P,csv_fh,row.names=FALSE)
            #qobj = qvalue(gwas_p_df$P,fdr.level = fdr_level)
            b = 0:20*.05
            png(filename = graph_fh)
            #hist(qobj)
            hist(gwas_p_df$P,breaks = b)
            dev.off()

            gwas_keep = qvalue(gwas_p_df$P,fdr.level = fdr_level)$significant
            """
            @rget(gwas_keep)
        catch
            println("Using Benjamini-Hochberg Method")
            @warn "Using Benjamini-Hochberg Method"

            @rput(gwas_p_df,fdr_level,gwas_fh)
            R"""
            library(qvalue)
            graph_fh = paste(gwas_fh,"_pval_graph.png",sep = "")
            csv_fh = paste(gwas_fh,"_pval.csv",sep = "")
            print(graph_fh)
            write.csv(gwas_p_df$P,csv_fh,row.names=FALSE)
            #qobj = qvalue(gwas_p_df$P,fdr.level = fdr_level)
            b = 0:20*.05
            png(filename = graph_fh)
            #hist(qobj)
            hist(gwas_p_df$P,breaks = b)
            dev.off()

            gwas_keep = qvalue(gwas_p_df$P,fdr.level = fdr_level,pi = 1)$significant
            """
            @rget(gwas_keep)
        end
    end

    if sum(gwas_keep) > 0
        gwas_keep_df = gwas_df[gwas_keep,:] # subset gwas file to only those that are significant
        CSV.write(gwas_keep_fh,gwas_keep_df) # write this file to disk
        return 1
    else
        println("No significant p-values at an FDR of $fdr_level")
        return 0
    end
end

function gwas_bonf(gwas_fh,gwas_keep_fh,sig)
    gwas_df = CSV.read(gwas_fh, header=false, DataFrame)
    rename!(gwas_df,["CHR","SNP","POS","A1","A2","N","AF1","BETA","SE","P"]) # assume this format of gwas result file
    gwas_p_df = select(gwas_df, :P)
    bonf = sig/size(gwas_p_df)[1]

    gwas_keep_df = gwas_df[gwas_df.P .< bonf,:]

    if size(gwas_keep_df)[1] > 0
        CSV.write(gwas_keep_fh,gwas_keep_df) # write this file to disk
        return 1
    else
        println("No Bonferroni significant p-values")
        return 0
    end
end

function gwas_extract_snps(gwas_fh,gwas_keep_fh,keep_snp_set,delim)
    # extract keep_snp_set of snps from a gwas file
    gwas_io = GZip.open(gwas_fh)
    gwas_keep_io = open(gwas_keep_fh,"w")
    i = 1
    for line in eachline(gwas_io)
        snp = split(line,delim)[2]
        if in(snp,keep_snp_set)
            write(gwas_keep_io,line*"\n")
        end
        i += 1
        if (i % 1000000) == 0
            #println(i)
        end
    end
    close(gwas_io)
    close(gwas_keep_io)
end

function gwas_innerjoin(gwas_fh_a,snp_col_a,gwas_fh_b,snp_col_b,out_fh)
    ### assume files have headers, wlll append names with _a and _b

    # load into DataFrames
    gwas_a_df = CSV.read(gwas_fh_a, DataFrame);
    gwas_b_df = CSV.read(gwas_fh_b, DataFrame);

    #rename cols
    rename!(gwas_a_df, names(gwas_a_df) .* "_a");
    rename!(gwas_a_df, names(gwas_a_df)[snp_col_a] => "SNP");

    rename!(gwas_b_df, names(gwas_b_df) .* "_b")
    rename!(gwas_b_df, names(gwas_b_df)[snp_col_b] => "SNP")

    gwas_joined_df = innerjoin(gwas_a_df,gwas_b_df, on = :SNP)

    CSV.write(out_fh,gwas_joined_df)
end

function gwas_shared_loci(somaid_a,somaid_b,fdr_level,root_path,met_flag)
    # determine shared loci between two gwas files

    data_folder = "$(root_path)/data/gwas_protein_out/"
    gwas_file_a_pre = "jhs_prot_adjust_age_sex_pcs_prot_"
    gwas_file_a_post = "_fdr_$(fdr_level)_cis.txt"
    gwas_file_b_pre = "jhs_prot_adjust_age_sex_pcs_prot_"
    gwas_file_b_post = "_keep.fastgwa.txt"

    if met_flag == 1

    end

    a_fh = data_folder*gwas_file_a_pre*somaid_a*gwas_file_a_post
    b_fh = data_folder*gwas_file_b_pre*somaid_b*gwas_file_b_post

    a_info = stat(a_fh)
    if a_info.size == 0
        return 0
    end
    b_info = stat(b_fh)
    if b_info.size == 0
        return 0
    end

    a_df = CSV.read(a_fh, DataFrame);
    a_snp_set = Set(a_df.SNP);

    b_sub_fh = "$(root_path)/data/gwas_protein_out/b_$(somaid_a)_$(somaid_b)_sub_temp.txt"

    gwas_extract_snps(b_fh,b_sub_fh,a_snp_set,"\t")

    b_keep_fh = "$(root_path)/data/gwas_protein_out/b_$(somaid_a)_$(somaid_b)_shared_temp.txt"

    any_shared_loci = gwas_fdr(b_sub_fh,b_keep_fh,fdr_level)

    if any_shared_loci == 1
        b_keep_df = CSV.read(b_keep_fh, DataFrame)
        a_sub_fh = "$(root_path)/data/gwas_protein_out/a_$(somaid_a)_$(somaid_a)_sub_temp.txt"
        shared_snp_set = Set(b_keep_df.SNP)
        gwas_extract_snps(a_fh,a_sub_fh,shared_snp_set,",")

        a_keep_df = CSV.read(a_sub_fh;header=false, DataFrame)
        DataFrames.rename!(a_keep_df, names(b_keep_df) .*"_a")
        DataFrames.rename!(b_keep_df, names(b_keep_df) .*"_b")
        shared_loci = hcat(a_keep_df,b_keep_df)

        shared_loci.SOMAID_a = somaid_a
        shared_loci.SOMAID_b = somaid_b
        shared_loci.BETA_b_BETA_a = shared_loci.BETA_b ./ shared_loci.BETA_a

        rm(a_sub_fh)
        rm(b_sub_fh)
        rm(b_keep_fh)

        return shared_loci
        #return 1
    else
        return 0
    end
end

function gwas_gen_ld_subdict(gwas_df,snp_col,ld_dict,ld_thresh)
    # generate a sub-dictionary of ld given a gwas dataframe and ld_threshold

    snp_arr = gwas_df[:,snp_col]
    gwas_ld_dict = Dict()
    for i in 1:length(snp_arr)
        for j in 1:length(snp_arr)
            if i != j
                temp_snp_set = Set([snp_arr[i],snp_arr[j]])
                if haskey(ld_dict,temp_snp_set) && ld_dict[temp_snp_set] > ld_thresh
                    gwas_ld_dict[temp_snp_set] = ld_dict[temp_snp_set]
                #else
                #    gwas_ld_dict[temp_snp_set] = 0
                end
            end
        end
    end

    return gwas_ld_dict
end

function gwas_ld_clump(gwas_fh,snp_col,p_col,ld_dict,ld_thresh)
    gwas_df = CSV.read(gwas_fh, DataFrame);
    sort!(gwas_df,p_col)
    
    gwas_ld_dict = gwas_gen_ld_subdict(gwas_df,snp_col,ld_dict,ld_thresh)

    println("Starting LD Clump...")
    for i in 1:size(gwas_df)[1]
        gwas_df_n = size(gwas_df)[1]
        if i >= gwas_df_n
            ## reached the end of the modified gwas_df
            break
        end

        for j in gwas_df_n:-1:(i+1)
            temp_snp_set = Set([gwas_df[i,snp_col],gwas_df[j,snp_col]])
            if haskey(ld_dict,temp_snp_set) && ld_dict[temp_snp_set] > ld_thresh
                delete!(gwas_df,j)
            end
        end
    end
    println("Finished!!!")

    return gwas_df
end

function gwas_ld_r(root_path,ldr_file,study)
    #remove "tss_tol" parameter

    # create list of all snps after ld clumping
    println("Starting gwas_ld_dict...")
    println("Create List of all snps after ld clumping...")
    cut_str = "cut -f2 -d, $(root_path)/data/gwas_protein_out/ldr/$(ldr_file) | grep \"chr\" > $(root_path)/data/gwas_protein_out/ldr/$(ldr_file).ldr_snps.txt"
    
    run(`bash -c "$(cut_str)"`)
    #run(cut_str)

    println("Finished!!!")

    # make bed file of only ldr_snps
    println("Make bed files of only ldr_snps...")
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

    plink_extract_str = `plink --bfile $(bfile_fh) --extract $(root_path)/data/gwas_protein_out/ldr/$(ldr_file).ldr_snps.txt --make-bed --keep-allele-order --out $(bfile_fh)_$(ldr_file).ldr_snps`
    run(plink_extract_str)

    println("Finished!!!")

    # compute ld file from new bed file
    #plink_ld_window = tss_tol / 1000
    #plink_ld_str = "plink --bfile $(bfile_fh)_$(ldr_file).ldr_snps --r --ld-window $(tss_tol) --ld-window-kb $(plink_ld_window) --ld-window-r2 0 --out $(root_path)/data/gwas_protein_out/ldr/$(ldr_file).ldr_snps_plink"

    plink_ld_str = "plink --bfile $(bfile_fh)_$(ldr_file).ldr_snps --r square yes-really --keep-allele-order --out $(root_path)/data/gwas_protein_out/ldr/$(ldr_file).ldr_snps_plink"
    run(`bash -c "$(plink_ld_str)"`)

end

function gwas_mr(gwasjoin_df,exp_name,out_name,use_corr,ld_r,psi) ### change to pass ld_r
    #perform mr using R package from an ld reduced gwas dataframe

    if use_corr == 1
        ## Do MR with correlation matrix

        ### Orient Beta's so that the effect always has a positive Beta
        orient_arr = []
        AF1_a_arr = gwasjoin_df.AF1_a
        for i in 1:length(AF1_a_arr)
            if sign(gwasjoin_df.BETA_a[i]) == -1
                gwasjoin_df.BETA_a[i] *= -1
                gwasjoin_df.BETA_b[i] *= -1
                push!(orient_arr,-1)
            else
                push!(orient_arr,1)
            end
        end

        orient_mat = orient_arr * transpose(orient_arr)
        ld_r = ld_r .* orient_mat
        
        @rput(gwasjoin_df,exp_name,out_name,ld_r)
        R"""

        require('MendelianRandomization')

        mr_input_obj = mr_input(gwasjoin_df$BETA_a,gwasjoin_df$SE_a,gwasjoin_df$BETA_b,gwasjoin_df$SE_b,exposure=exp_name,outcome=out_name, corr = as.matrix(ld_r))

        if(length(mr_input_obj$snps) >= 3000){

            mrv = mr_allmethods(mr_input_obj, method = "all")$Values
            
            mrv$name = paste(exp_name,"->",out_name,sep = "")

            names(mrv) = c("Method","B","SE","95CI_L","95CI_U","P","name")

            # reorder rows to put IVW first
            mrv = mrv[c(4:7,1:3,8:15),]

            mrv$Method[c(9,11,13,15)]

            mrv$Method[c(9,11,13,15)] = c("MRE_intercept","P_MRE_intercept","R_MRE_intercept","PR_MRE_intercept")

        }else {
            mrv = mr_allmethods(mr_input_obj, method = "ivw")$Values
            
            mrv$name = paste(exp_name,"->",out_name,sep = "")
        
            names(mrv) = c("Method","B","SE","95CI_L","95CI_U","P","name")
        }

        mrv = reshape(mrv,idvar = "name", timevar = "Method", direction="wide")
        
        """
    else
        println("DO MR w/out correlation matrix")
        ## DO MR without correlation matrix
        @rput(gwasjoin_df,exp_name,out_name,psi)
        R"""

        require('MendelianRandomization')
        print("Create MR Input Object")
        mr_input_obj = mr_input(gwasjoin_df$BETA_a,gwasjoin_df$SE_a,gwasjoin_df$BETA_b,gwasjoin_df$SE_b,exposure=exp_name,outcome=out_name)

        if(length(mr_input_obj$snps) >= 3){

            #mrv = mr_allmethods(mr_input_obj, method = "all")$Values
            print("Start main MR")
            mrv = mr_allmethods(mr_input_obj, method = "main")$Values
            print("Start MaxLik MR")
            if(psi != 999){
                mrv_liml = mr_maxlik(mr_input_obj,psi = psi)
                mrv = rbind(mrv,c("MaxLik",mrv_liml$Estimate,mrv_liml$StdError,mrv_liml$CILower,mrv_liml$CIUpper,mrv_liml$Pvalue))
            }

            #mrv$Method = c("med","med_wt","med_pen_wt","ivw","ivw_pen","ivw_rob","ivw_pen_rob","egg","egg_int","egg_pen","egg_pen_int","egg_rob","egg_rob_int","egg_pen_rob","egg_pen_rob_int")
            
            ### test dimensions of gwasjoin_df to determine which methods should be in mrv$Method
            # can MR Median be done with single one??? probably

            ##mrv$Method = c("med","med_wt","ivw","egg","egg_int")
            
            mrv$name = paste(exp_name,"->",out_name,sep = "")

            names(mrv) = c("Method","B","SE","95CI_L","95CI_U","P","name")
            mrv$Method[5] = "MRE_intercept"

            # reorder rows to put IVW first
            #mrv = mrv[c(4:7,1:3,8:15),]

            #mrv$Method[c(9,11,13,15)]

            #mrv$Method[c(9,11,13,15)] = c("MRE_intercept","P_MRE_intercept","R_MRE_intercept","PR_MRE_intercept")


            # reorder rows to put IVS and MaxLik first
            mrv = mrv[c(3,6,1,2,4,5),]

        }else {
            #mrv = mr_allmethods(mr_input_obj, method = "ivw")$Values
            print("Start IVW MR")
            mr_ivw = mr_ivw(mr_input_obj)
            mrv = as.data.frame(t(c("ivw",mr_ivw$Estimate,mr_ivw$StdError,mr_ivw$CILower,mr_ivw$CIUpper,mr_ivw$Pvalue)))
            print("Start MaxLik MR")
            if(psi != 999){
                mrv_liml = mr_maxlik(mr_input_obj,psi = psi)
                mrv = rbind(mrv,c("MaxLik",mrv_liml$Estimate,mrv_liml$StdError,mrv_liml$CILower,mrv_liml$CIUpper,mrv_liml$Pvalue))
            }

            
            mrv$name = paste(exp_name,"->",out_name,sep = "")

            names(mrv) = c("Method","B","SE","95CI_L","95CI_U","P","name")
        }

        mrv = reshape(mrv,idvar = "name", timevar = "Method", direction="wide")
        
        """
    end

    println("finished mr function")

    return @rget(mrv)
end
