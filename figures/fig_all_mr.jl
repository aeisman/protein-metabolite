using CSV,DataFrames,RCall,StatsBase
try
    include("../analysis/omics_tools.jl")
    include("../analysis/gwas_tools.jl")
catch
    include("./analysis/omics_tools.jl")
    include("./analysis/gwas_tools.jl")
end

### Load config ###
conf = "./figures/fig.conf";
project_path = read_conf(conf,"project_path");
results_sub_path = read_conf(conf,"results_sub_path");
mr_data_file = read_conf(conf,"mr_data_file");

cor_data_file = read_conf(conf,"cor_data_file");

cor_sig_enriched_union_file = read_conf(conf,"cor_sig_enriched_union_file");
master_snp_maf_union_gwas_file = read_conf(conf,"master_snp_maf_union_gwas_file")
gene_locs_file = read_conf(conf,"gene_locs_file")

fig_mr_out = read_conf(conf,"fig_mr_out")

annot_sub_path = read_conf(conf,"annot_sub_path")
met_anno_fh = project_path*annot_sub_path*"ogm_meta_key_anno.csv"

#    met_anno = CSV.read(met_anno_fh, DataFrame)
    

### Load SNP data ###

### Load correlation data ###
cor_df = CSV.read(project_path*results_sub_path*cor_data_file,DataFrame)

### Load MR results ###
mr = CSV.read(project_path*results_sub_path*mr_data_file,DataFrame)
mr[!,:sig] .= 0;
mr[mr.META_P_ivw .< 0.05,:sig] .= 1;
mr[mr.META_ivw_qvalue .< 0.1,:sig] .= 2;
mr[mr.META_ivw_qvalue .< 0.05,:sig] .= 3;
sort(countmap(mr.sig))

mr_sig = filter(row -> row.sig == 3, mr)

### Load correlation and enrichment results
cor_sig_enriched_union_df = CSV.read(project_path*results_sub_path*cor_sig_enriched_union_file,DataFrame)
cor_sig_enriched_union_df[!,:soma_id] .= ""
for i in 1:size(cor_sig_enriched_union_df)[1]
    cor_sig_enriched_union_df.soma_id[i] = match(r"SL[0-9]+",cor_sig_enriched_union_df.cor_id[i]).match
end

### Create figures

# Proteins w/ CIS instruments
n_cis = length(Set(mr[:,Symbol("p.soma_id")]))
n_all = length(Set(cor_sig_enriched_union_df.soma_id))
n_no_cis = n_all - n_cis

bar_table = [n_cis,n_no_cis]
n_instruments_fh = project_path*fig_mr_out*"n_instruments.pdf"

@rput(bar_table,n_all,n_instruments_fh)
R"""
    pdf(n_instruments_fh, width = 3, height = 7)
    par(mar = c(7,5,0,3))
    barplot(as.matrix(bar_table), beside = FALSE, ylim = c(0,1400), cex.axis = 2, cex.lab = 2, ylab = "", col = c("dodgerblue2","white"), axes = FALSE)

    axis(2,at=c(0,bar_table[1],n_all),cex.axis = 2)
    #text(0.7,1365,n_all, cex = 2)
    dev.off()
"""

# By Chromosome
snps = CSV.read(project_path*results_sub_path*master_snp_maf_union_gwas_file,DataFrame)
snps_select = filter(row -> row.jhs_select == 1 || row.mesa_select == 1 || row.her_select == 1, snps)

gene_locs = CSV.read(project_path*annot_sub_path*gene_locs_file,DataFrame)

snps_select[!,:chr] .= ""
for i in 1:size(snps_select)[1]
    t_chr = match(r"(?<=chr)\d?\d(?=:)",snps_select.snp[i]).match
    snps_select.chr[i] = t_chr
end
chr_ct_dict = countmap(snps_select.chr)
chr_ct_arr = []
for i in 1:22
    push!(chr_ct_arr,chr_ct_dict["$(i)"])
end

mr_prots = Set(mr[:,Symbol("p.soma_id")])

chr_p_ct_dict = countmap(parse.(Int,filter(row -> in(row.somaid,mr_prots), gene_locs).chromosome))
chr_p_ct_arr = []
for i in 1:22
    push!(chr_p_ct_arr,chr_p_ct_dict[i])
end

chr_ct_fh = project_path*fig_mr_out*"chr_ct.pdf"
chr_p_ct_fh = project_path*fig_mr_out*"chr_p_ct.pdf"

@rput(chr_ct_arr,chr_p_ct_arr,chr_ct_fh,chr_p_ct_fh)
R"""
    pdf(chr_ct_fh, width = 4, height = 7)
    
    names_arr = 1:22
    barplot(t(chr_ct_arr),horiz = TRUE, ylab = "Chromosome", xlab = "MR Instruments (n)", names.arg = names_arr, xlim = c(0,200), cex.axis = 1.25, cex.lab = 1.5, col = c("dodgerblue2"))
    
    dev.off()

    
    pdf(chr_p_ct_fh, width = 4, height = 7)
    
    names_arr = 1:22
    barplot(t(chr_p_ct_arr),horiz = TRUE, ylab = "Chromosome", xlab = "Proteins (n)", names.arg = names_arr, xlim = c(0,75), cex.axis = 1.25, cex.lab = 1.5, col = c("dodgerblue2"))
    
    dev.off()
"""

# SNPs proximity to test_scatter

snps_select[!,:pos] .= ""
for i in 1:size(snps_select)[1]
    snps_select.pos[i] = split(snps_select.snp[i],":")[2]
end

tss_dict = Dict()
for i in 1:size(gene_locs)[1]
    tss_dict[gene_locs.somaid[i]] = gene_locs.tss[i]
end

snps_select[!,:pos_m_tss] .= 999999999
for i in 1:size(snps_select)[1]
    try
        snps_select.pos_m_tss[i] = parse(Int64,snps_select.pos[i]) - tss_dict[snps_select.prot[i]]
    catch
        @warn("No Gene Locs Info for $(snps_select.prot[i])")
    end
end

pos_m_tss_arr = snps_select.pos_m_tss
l = 10^6
for i in 1:length(pos_m_tss_arr)
    if pos_m_tss_arr[i] > l
        pos_m_tss_arr[i] = l
    elseif pos_m_tss_arr[i] < -l
        pos_m_tss_arr[i] = -l
    end
end

pos_m_tss_fh = project_path*fig_mr_out*"pos_m_tss.pdf"
pos_m_tss_density_fh = project_path*fig_mr_out*"pos_m_tss_density.pdf"

@rput(pos_m_tss_arr,l,pos_m_tss_fh,pos_m_tss_density_fh)
R"""
    pdf(pos_m_tss_fh, width = 4, height = 6)
    hist(pos_m_tss_arr, xlim = c(-l,l), ylim = c(0,600), breaks = 10)
    dev.off()
    
    pdf(pos_m_tss_density_fh, width = 4, height = 6)
    par(mar = c(5,0,0,0))
    plot(density(pos_m_tss_arr), lwd = 5, cex.lab = 2, xlab = "", ylab = "", main = "", axes = FALSE, col = c("dodgerblue2"))
    axis(1, at = c(-1000000,-500000,0,500000,1000000), labels = c("-1M","-500K","0","500K","1M"))
    dev.off()

"""

# Volcano Plot
function volcano_plot(df,p_col,b_col,sig_var,fn)
    #sort!(df,Symbol(sig_var))
    df = filter(row -> !ismissing(row[Symbol(p_col)]), df)
    
    log10_p = -log.(10,df[:,Symbol(p_col)])
    betas = df[:,Symbol(b_col)]
    #sig = (df[:,Symbol(sig_var)] .== 1) .+ 1
    sig_arr = df[:,Symbol(sig_var)]
    color = []
    shape = []
    for i in 1:length(sig_arr)
        if df[:,Symbol(p_col)][i] > 0.05
            push!(color,"gray50")
            push!(shape,"1")
        elseif sig_arr[i] == 1
            push!(color,"green3")
            push!(shape,"16")
        elseif sig_arr[i] == 0
            push!(color,"black")
            push!(shape,"1")
        end
    end
    v_plot_data = hcat(betas,log10_p)
    v_plot_data = hcat(v_plot_data,color)
    v_plot_data = hcat(v_plot_data,shape)

    v_plot_data = DataFrame(betas = betas, log10_p = log10_p, color = color, shape = shape)

    v_plot_data_not_sig = filter(row -> row.color == "gray50")
    v_plot_data_
    @rput(v_plot_data,fn)
    R"""
        pdf(fn)

        plot(v_plot_data$log10_p~v_plot_data$betas,xlim = c(-max(abs(v_plot_data$betas)),max(abs(v_plot_data$betas))), pch = 6, xlab = "META-Analyzed Beta Coefficient", ylab = "-log10(q-value)")
        abline(h = 1.3, lty = 2, lwd = 2, col = "gray50")
        #points(v_plot_data$log10_p~v_plot_data$betas, col = v_plot_data$color, pch = v_plot_data$shape)
        

        dev.off() 
    """
end

#df = sort(apoe_cor_df,:var2_Combined_Lipids,rev=true)
#p_col = "meta_p_qval"
#b_col = "meta_rho"
#sig_var = "var2_Combined_Lipids"
#fn = "apoe_lipid_volcano.pdf"
#volcano_plot(df, p_col, b_col, sig_var, fn)

function volcano_plot_split(df,p_col,b_col,sig_var,fn,y_split)
    #y_split = 10
    #df = filter(row -> row.META_sig_qvalue < 0.07,mr_filt_sig)
    #df = mr_filt_sig
    #p_col = "META_sig_qvalue"
    #b_col = "META_B_ivw"
    #sig_var = "META_sig"
    #fn = "mr_filt_split_volcano.pdf"

    sort!(df,Symbol(sig_var))
    
    log10_p = -log.(10,df[:,Symbol(p_col)])
    betas = df[:,Symbol(b_col)]
    sig_var_col = df[:,Symbol(sig_var)]
    if typeof(sig_var_col[1]) == String
        #for i in 1:length(sig_var_col)
        #    if sig_var_col[i] == "grey"
        #        sig[i] = "Gray50"
        #    elseif sig_var_col[i] == "black"
        #        sig[i] = "black"
        #    elseif sig_var_col[i] == "red"
        #        sig[i] == "red"
        #    end
        #end
        sig = sig_var_col
    else
        sig = (df[:,Symbol(sig_var)] .== 3) .+ 1
    end
    
    v_plot_data = hcat(betas,log10_p)
    v_plot_data = hcat(v_plot_data,sig)

    v_plot_data = DataFrame(beta = betas, log10_p = log10_p, sig = sig)

    y_max = maximum(log10_p)

    v_plot_data_upper = filter(row -> row.log10_p > 10, v_plot_data)
    v_plot_data_lower = filter(row -> row.log10_p <= 10, v_plot_data)

    @rput(v_plot_data_lower,v_plot_data_upper,v_plot_data,fn,y_max,y_split)
    R"""
        cnvrt.coords <-function(x,y=NULL){
        # Stolen from the teachingDemos library, simplified for this use case
            xy <- xy.coords(x,y, recycle=TRUE)
            cusr <- par('usr')
            cplt <- par('plt')	
            plt <- list()
            plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
            plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])
            fig <- list()
            fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
            fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
            return( list(fig=fig) )
        }
        subplot <- function(fun, x, y=NULL){
            # Stolen from the teachingDemos library, simplified for this use case
            old.par <- par(no.readonly=TRUE)
            on.exit(par(old.par))
            xy <- xy.coords(x,y)
            xy <- cnvrt.coords(xy)$fig
            par(plt=c(xy$x,xy$y), new=TRUE)
            fun
            tmp.par <- par(no.readonly=TRUE)
            return(invisible(tmp.par))
        }

        #fn = "mr_volcano.pdf"
        pdf(fn, width = 8, height = 5)
        par(mar=c(3, 6, 3, 3))
        x_limits = c(-max(abs(v_plot_data[,1])),max(abs(v_plot_data[,1])))
        x_limits = c(-1.5,1.5)
        y_lower = c(0,y_split)
        y_upper = c(y_split,max(v_plot_data$log10_p))
        y_upper = c(y_split,40)

        plot(c(0,1), c(0,11), type='n', axes=FALSE, ylab="-log10(q-value)", xlab='Beta Coefficient')
        y_sig = 1.3
        segments(x0 = 0.05, x1 = 0.95, y0 = y_sig, y1 = y_sig, lty = 2, lwd = 2, col = "gray50")

        subplot(plot(v_plot_data_lower[,2]~v_plot_data_lower[,1], xlim = x_limits, ylim = y_lower, col = v_plot_data_lower[,3], xlab = "", ylab = "", frame = FALSE), x=c(0,1), y=c(0,8))

        subplot(plot(v_plot_data_upper[,2]~v_plot_data_upper[,1], xlim = x_limits, ylim = y_upper, col = v_plot_data_upper[,3], xlab = "", ylab = "", xpd=FALSE, xaxt = "n", frame = FALSE), x=c(0,1), y=c(8.25,10))

        dev.off() 
    """
end

function volcano_plot_split2(df,p_col,b_col,sig_var,fn,y_split)
    #y_split = 10
    #df = filter(row -> row.META_sig_qvalue < 0.07,mr_filt_sig)
    #df = mr_filt_sig
    #p_col = "META_sig_qvalue"
    #b_col = "META_B_ivw"
    #sig_var = "META_sig"
    #fn = "mr_filt_split_volcano.pdf"

    sort!(df,Symbol(sig_var))

    df[!,:log10_p] = -log.(10,df[:,Symbol(p_col)])

    mouse_prots = Set(["SL000276","SL005574","SL000668"])

    sig_df = filter(row -> !in(row[Symbol("p.soma_id")],mouse_prots) && row[Symbol(sig_var)] == 3,df)
    sig_df_lower = filter(row -> row.log10_p < 10, sig_df)
    sig_df_upper = filter(row -> row.log10_p > 10, sig_df)

    not_sig_df = filter(row -> row[Symbol(sig_var)] != 3,df)

    # ApoE SL000276
    apoe_df = filter(row -> row[Symbol("p.soma_id")] == "SL000276" && row[Symbol(sig_var)] == 3,df)
    apoe_df_lower = filter(row -> row.log10_p < 10, apoe_df)
    apoe_df_upper = filter(row -> row.log10_p > 10, apoe_df)


    # Acy-1 SL005574
    acy1_df = filter(row -> row[Symbol("p.soma_id")] == "SL005574" && row[Symbol(sig_var)] == 3, df)
    acy1_df_lower = filter(row -> row.log10_p < 10, acy1_df)
    acy1_df_upper = filter(row -> row.log10_p > 10, acy1_df)
    

    # CD36 SL000668
    cd36_df = filter(row -> row[Symbol("p.soma_id")] == "SL000668" && row[Symbol(sig_var)] == 3, df)
    cd36_df_lower = filter(row -> row.log10_p < 10, cd36_df)
    cd36_df_upper = filter(row -> row.log10_p > 10, cd36_df)
    
    
    log10_p = -log.(10,df[:,Symbol(p_col)])
    betas = df[:,Symbol(b_col)]
    sig_var_col = df[:,Symbol(sig_var)]
    if typeof(sig_var_col[1]) == String
        #for i in 1:length(sig_var_col)
        #    if sig_var_col[i] == "grey"
        #        sig[i] = "Gray50"
        #    elseif sig_var_col[i] == "black"
        #        sig[i] = "black"
        #    elseif sig_var_col[i] == "red"
        #        sig[i] == "red"
        #    end
        #end
        sig = sig_var_col
    else
        sig = (df[:,Symbol(sig_var)] .== 3) .+ 1
    end
    
    v_plot_data = hcat(betas,log10_p)
    v_plot_data = hcat(v_plot_data,sig)

    v_plot_data = DataFrame(beta = betas, log10_p = log10_p, sig = sig)

    y_max = maximum(log10_p)

    v_plot_data_upper = filter(row -> row.log10_p > 10, v_plot_data)
    v_plot_data_lower = filter(row -> row.log10_p <= 10, v_plot_data)

    #@rput(v_plot_data_lower,v_plot_data_upper,v_plot_data,fn,y_max,y_split)
    @rput(v_plot_data,sig_df_lower, sig_df_upper, apoe_df_lower, apoe_df_upper, acy1_df_lower, acy1_df_upper, cd36_df_lower, cd36_df_upper, not_sig_df, fn, y_max, y_split)
    R"""
        cnvrt.coords <-function(x,y=NULL){
        # Stolen from the teachingDemos library, simplified for this use case
            xy <- xy.coords(x,y, recycle=TRUE)
            cusr <- par('usr')
            cplt <- par('plt')	
            plt <- list()
            plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
            plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])
            fig <- list()
            fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
            fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
            return( list(fig=fig) )
        }
        subplot <- function(fun, x, y=NULL){
            # Stolen from the teachingDemos library, simplified for this use case
            old.par <- par(no.readonly=TRUE)
            on.exit(par(old.par))
            xy <- xy.coords(x,y)
            xy <- cnvrt.coords(xy)$fig
            par(plt=c(xy$x,xy$y), new=TRUE)
            fun
            tmp.par <- par(no.readonly=TRUE)
            return(invisible(tmp.par))
        }

        #fn = "mr_volcano.pdf"
        pdf(fn, width = 8, height = 5)
        par(mar=c(3, 6, 3, 3))
        x_limits = c(-max(abs(v_plot_data[,1])),max(abs(v_plot_data[,1])))
        x_limits = c(-1.5,1.5)
        y_lower = c(0,y_split)
        y_upper = c(y_split,max(v_plot_data$log10_p))
        y_upper = c(y_split,40)

        plot(c(0,1), c(0,11), type='n', axes=FALSE, ylab="-log10(q-value)", xlab='Beta Coefficient')
        y_sig = 1.3
        segments(x0 = 0.05, x1 = 0.95, y0 = y_sig, y1 = y_sig, lty = 2, lwd = 2, col = "gray50")

        subplot(plot(not_sig_df$log10_p~not_sig_df$META_B_ivw, xlim = x_limits, ylim = y_lower, xlab = "", ylab = "", frame = FALSE), x=c(0,1), y=c(0,8))
        subplot(plot(apoe_df_lower$log10_p~apoe_df_lower$META_B_ivw, xlim = x_limits, ylim = y_lower, xlab = "", ylab = "", frame = FALSE, add = TRUE, pch = 15, col = "red"), x=c(0,1), y=c(0,8))
        subplot(plot(acy1_df_lower$log10_p~acy1_df_lower$META_B_ivw, xlim = x_limits, ylim = y_lower, xlab = "", ylab = "", frame = FALSE, add = TRUE, pch = 17, col = "red"), x=c(0,1), y=c(0,8))
        subplot(plot(cd36_df_lower$log10_p~cd36_df_lower$META_B_ivw, xlim = x_limits, ylim = y_lower, xlab = "", ylab = "", frame = FALSE, add = TRUE, pch = 18, col = "red"), x=c(0,1), y=c(0,8))
        subplot(plot(sig_df_lower$log10_p~sig_df_lower$META_B_ivw, xlim = x_limits, ylim = y_lower, xlab = "", ylab = "", frame = FALSE, add = TRUE, col = "red"), x=c(0,1), y=c(0,8))

        subplot(plot(apoe_df_upper$log10_p~apoe_df_upper$META_B_ivw, xlim = x_limits, ylim = y_upper, xlab = "", ylab = "", xpd=FALSE, xaxt = "n", frame = FALSE, add = TRUE, pch = 15, col = "red"), x=c(0,1), y=c(8.25,10))
        subplot(plot(acy1_df_upper$log10_p~acy1_df_upper$META_B_ivw, xlim = x_limits, ylim = y_upper, xlab = "", ylab = "", xpd=FALSE, xaxt = "n", frame = FALSE, add = TRUE, pch = 17, col = "red"), x=c(0,1), y=c(8.25,10))
        subplot(plot(cd36_df_upper$log10_p~cd36_df_upper$META_B_ivw, xlim = x_limits, ylim = y_upper, xlab = "", ylab = "", xpd=FALSE, xaxt = "n", frame = FALSE, add = TRUE, pch = 18, col = "red"), x=c(0,1), y=c(8.25,10))
        subplot(plot(sig_df_upper$log10_p~sig_df_upper$META_B_ivw, xlim = x_limits, ylim = y_upper, xlab = "", ylab = "", xpd=FALSE, xaxt = "n", frame = FALSE, add = TRUE, col = "red"), x=c(0,1), y=c(8.25,10))
#, col = "red")
        #subplot(plot(v_plot_data_upper[,2]~v_plot_data_upper[,1], xlim = x_limits, ylim = y_upper, col = v_plot_data_upper[,3], xlab = "", ylab = "", xpd=FALSE, xaxt = "n", frame = FALSE), x=c(0,1), y=c(8.25,10))
        #subplot(plot(v_plot_data_upper[,2]~v_plot_data_upper[,1], xlim = x_limits, ylim = y_upper, col = v_plot_data_upper[,3], xlab = "", ylab = "", xpd=FALSE, xaxt = "n", frame = FALSE), x=c(0,1), y=c(8.25,10))
        #subplot(plot(v_plot_data_upper[,2]~v_plot_data_upper[,1], xlim = x_limits, ylim = y_upper, col = v_plot_data_upper[,3], xlab = "", ylab = "", xpd=FALSE, xaxt = "n", frame = FALSE), x=c(0,1), y=c(8.25,10))
        #subplot(plot(v_plot_data_upper[,2]~v_plot_data_upper[,1], xlim = x_limits, ylim = y_upper, col = v_plot_data_upper[,3], xlab = "", ylab = "", xpd=FALSE, xaxt = "n", frame = FALSE), x=c(0,1), y=c(8.25,10))

        legend(x = "topright", legend = c("Acy1 (6 metabolites)", "ApoE (10 metabolites)", "CD36 (53 metabolites)", "Other (49 proteins, 90 metabolites)","Not Significant"),pch = c(17,15,18,1,1), col = c("red","red","red","red","black"), bty = "n")


        dev.off() 
    """
end

#p_col = "META_P_ivw"
p_col = "META_ivw_qvalue"
b_col = "META_B_ivw"
sig_var = "sig"
y_split = 10
#fh = "/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_figures/pm_figures_current_indesign_files/fig_all_mr_files/mr_volcano_2.pdf"
fh = project_path*fig_mr_out*"mr_volcano_2.pdf"
volcano_plot_split2(mr,p_col,b_col,sig_var,fh,y_split)
volcano_plot_split(mr_filt_sig,p_col,b_col,sig_var,"mr_volcano.pdf",y_split)


# load mouse data
mouse = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/data/mr_filt_sig_final_20220202 for mouse coverage volcano.csv",DataFrame)
sig_var = "color"
volcano_plot_split(mouse,p_col,b_col,sig_var,"mr_volcano_mouse.pdf",y_split)




## WHAT CHANGED?

mr_filt_sig_05 = filter(row -> row.META_sig == 3,mr_filt_sig)
mr_old = CSV.read("/Users/aaroneisman/Downloads/v4_mr_pm_qval_by_protein_qval_0.05.csv",DataFrame)

old_set = Set(mr_old.cor_id)
new_set = Set(mr_filt_sig_05.cor_id)
in_common_set = intersect(old_set,new_set)

lost_set = setdiff(old_set,new_set)

filter(row -> in(row.cor_id,lost_set),mr_old)

lost_df = sort(filter(row -> in(row.cor_id,lost_set),mr_old),:p_name)

apoe_iso_set = Set(["SL000277","SL004668","SL004669"])
apoe_iso_mets = Set(filter(row -> in(row.p_soma_id,apoe_iso_set),lost_df).m_id)

apo_mr_df = filter(row -> row[Symbol("p.soma_id")] == "SL000276",mr_filt_sig_05)
apo_mets = Set(filter(row -> in(row.p_soma_id,apoe_iso_set),lost_df).m_id)

cd36_lost_set = Set(["OGM310793","OGM315409","OGM315226","OGM310014","OGM319851"])