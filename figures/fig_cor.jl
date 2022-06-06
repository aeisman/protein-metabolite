using DataFrames,CSV,RCall,StatsBase,Printf,Random
try
    include("../analysis/omics_tools.jl")
catch
    include("./analysis/omics_tools.jl")
end

function remove_dups(arr)
    # remove duplicates in an arr but maintain ordering
    new_arr_set = Set()
    new_arr = []

    for i in 1:length(arr)
        if !in(arr[i], new_arr_set)
            push!(new_arr,arr[i])
            push!(new_arr_set,arr[i])
        end
    end
    return new_arr
end

# Functions to process cor_df output from combine_correlations.jl

function clean_met_anno!(cor_df,met_anno)
    cor_df_mets = Set(cor_df.OGM)
    met_anno_mets = Set(met_anno.ogm_id)

    met_set = intersect(cor_df_mets,met_anno_mets)

    filter!(row -> in(row.ogm_id,met_set), met_anno)
end

function sort_by(df,order_var,order_df)
    order_arr = []

    for j in 1:size(order_df)[1]
        for i in 1:size(df)[1]
            if order_df[j,Symbol(order_var)] == df[i,Symbol(order_var)]
                push!(order_arr,i)
            end
        end
    end

    return order_arr
end

function gen_cor_mat(cor_df,met_anno)

    prots = sort(collect(Set(cor_df.PROT)))
    n_prots = length(prots)
    mets = met_anno.ogm_id

    n_mets = length(mets)

    cor_mat = zeros(n_prots,n_mets);
    cor_mat_sig = zeros(n_prots,n_mets);
    p_mat = ones(n_prots,n_mets)
    q_mat = ones(n_prots,n_mets)

    for i in 1:n_prots
        println(i)
        t1_cor_df = filter(row -> row.PROT == prots[i], cor_df)
        for j in 1:n_mets
            try
                t2_cor_df = filter(row -> row.OGM == met_anno.ogm_id[j], t1_cor_df)
                #t2_cor_df = filter(row -> row.OGM == mets[j], t1_cor_df)

                if size(t2_cor_df)[1] > 1
                    @warn("$(i) $(j) $(prots[i]) vs $(mets[j]) is not unique!!!")
                end

                r = t2_cor_df.meta_rho[1]
                p = t2_cor_df.meta_p[1]
                q = t2_cor_df.meta_p_qval[1]
                if r >= -1 && r <= 1
                    cor_mat[i,j] = r
                    p_mat[i,j] = p
                    q_mat[i,j] = q
                    if t2_cor_df.meta_p_qval[1] < 0.05
                        cor_mat_sig[i,j] = r
                    end
                end
            catch
                @warn("$(j) $(prots[i]) vs $(mets[j]) does not have a correlation value")
            end
        end
    end

    znorm_cor_mat = copy(cor_mat)
    for j in 1:size(cor_mat)[2]
        mean_abs_rho = mean(abs.(cor_mat[:,j]))
        std_rho = std(cor_mat[:,j])
        #cor_mat[:,j] = cor_mat[:,j] ./ mean_abs_rho
        znorm_cor_mat[:,j] = cor_mat[:,j] ./ std_rho
    end

    return cor_mat,znorm_cor_mat,cor_mat_sig,p_mat,q_mat,prots,mets
end

function trim_cor_mat(cor_mat,thresh)
    t_cor_mat = cor_mat[:,:];

    for i in 1:size(t_cor_mat)[1]
        for j in 1:size(t_cor_mat)[2]
            if t_cor_mat[i,j] < -thresh
                t_cor_mat[i,j] = -(thresh+0.01)
            elseif t_cor_mat[i,j] > thresh
                t_cor_mat[i,j] = (thresh+0.01)
            end
        end
    end

    return t_cor_mat
end

function find_hmap_point(hmap,prots,mets,prot,met)


    return x,y
end

### Functions to generate figures

function gen_cor_sig_bar_col(df,sig_var,rho_var,sig)
    t_df = DataFrame(p_qval = df[!,Symbol(sig_var)], rho = df[!,Symbol(rho_var)])
    dropmissing!(t_df)
    n_sig_neg = sum((t_df.p_qval .< sig) .* (t_df.rho .< 0))
    n_sig_pos = sum((t_df.p_qval .< sig) .* (t_df.rho .> 0))
    n_not_sig = length(t_df.rho) - n_sig_neg - n_sig_pos
    cor_sig_bar_col = [n_sig_neg,n_not_sig,n_sig_pos]
    return cor_sig_bar_col
end

function create_bargraphs(stacked_bar_table, out_fh)
    @rput(stacked_bar_table,out_fh)
    R"""
        pdf(out_fh, width = 10, height = 10)
        par(mar=c(14, 6, 3, 3))
        
        barplot(stacked_bar_table, col = c("navy","gray","firebrick3"), names.arg = c("JHS","MESA","HERITAGE","META"), ylab = expression(paste("Correlation Tests (n x ",10^{3},")")), ylim = c(0,500), cex.names = 2, cex.lab = 2, cex.axis = 2, legend.text = c("Negative Correlation", "No Correlation (q > 0.05)", "Positive Correlation"), args.legend = list(x = "bottom", inset = c(0,-.35), cex = 2))
    
        dev.off()
    """    
end

function create_adj_bargraphs(adj_bar_table, bmi_pct_str, egfr_pct_str, out_fh)
    @rput(adj_bar_table,bmi_pct_str,egfr_pct_str,out_fh)
    R"""
        pdf(out_fh, width = 10, height = 10)
        par(mar=c(14, 6, 3, 3))
        
        barplot(adj_bar_table, density = c(5,10,20), names.arg = c("Age/Sex","+BMI","+eGFR"), ylab = expression(paste("Significant Correlations (n x ",10^{3},")")), ylim = c(0,200), cex.names = 2, cex.lab = 2, cex.axis = 2)
    
        text(x = c(1.95,3.15), y = c(adj_bar_table[2]+10,adj_bar_table[3]+10), labels = c(bmi_pct_str,egfr_pct_str), cex = 2)
        
        dev.off()
    """    
end

function gen_scatter_df(df,var1,var2)
    t_df = df[:,[Symbol(var1),Symbol(var2)]]
    dropmissing!(t_df)
    return t_df
end

function create_scatter_plot(df,xname,yname,out_fh)
    #df = jhs_mesa_scatter_df
    #fn = "test_scatter.pdf"
    #xname = "JHS"
    #yname = "MESA"
    @rput(df,out_fh,xname,yname)
    R"""
    pdf(out_fh, width = 10, height = 10)
    par(mar=c(6, 6, 6, 6))

    plot(df, xlab = xname, ylab = yname, cex.lab = 2, cex.axis = 2, col = "gray50", xlim = c(-1,1), ylim = c(-1,1))
    abline(h = 0, lwd = 3)
    abline(v = 0, lwd = 3)
    abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)

    dev.off()
    """
end

function create_heatmap(t_cor_mat,met_anno,out_fh,rowv,colv)
    @rput(t_cor_mat, met_anno, out_fh, rowv, colv);
    R"""
        require('heatmap3')

        if (rowv == "NULL"){
            rowv = NULL
        }else if (rowv == "NA"){
            rowv = NA
        }

        if (colv == "NULL"){
            colv = NULL
        }else if (colv == "NA"){
            colv = NA
        }
    
        pdf(out_fh, width = 10, height = 10)    

        t_cor_mat = t(t_cor_mat)
        rcolors = cbind(pm_subclass = met_anno$pm_subclass_color,pm_spacer = "white", pm_class = met_anno$pm_class_color)
        hmap = heatmap3(t_cor_mat, method = "complete", Rowv = rowv, Colv = colv, balanceColor = TRUE, scale = "none", labRow = FALSE, labCol = FALSE, showColDendro = TRUE, RowSideColors = rcolors, col = colorRampPalette(c("navy", "white", "firebrick3"))(15))

        dev.off()
    """

    @rget(hmap)

    return(hmap)

end



function redo_bw_colors(met_anno)

    bw_colors = ["Black","Gray50"]

    met_anno[!,:pm_class_color] .= "white"
    met_anno[!,:pm_subclass_color] .= "white"

    pm_class = met_anno.pm_class[1]
    pm_subclass = met_anno.pm_subclass[1]
    t_pm_subclass_bin = false
    t_pm_class_bin = false
    for i in 1:size(met_anno)[1]
        t_pm_class = met_anno.pm_class[i]
        t_pm_subclass = met_anno.pm_subclass[i]

        if t_pm_class == pm_class
            met_anno.pm_class_color[i] = bw_colors[t_pm_class_bin + 1]
        else
            pm_class = t_pm_class
            t_pm_class_bin = !t_pm_class_bin
            met_anno.pm_class_color[i] = bw_colors[t_pm_class_bin + 1]
        end

        if t_pm_subclass == pm_subclass
            met_anno.pm_subclass_color[i] = bw_colors[t_pm_subclass_bin + 1]
        else
            pm_subclass = t_pm_subclass
            t_pm_subclass_bin = !t_pm_subclass_bin
            met_anno.pm_subclass_color[i] = bw_colors[t_pm_subclass_bin + 1]
        end

    end

    return met_anno
end

function do_fig_cor()

    # load and prep data and annotations

    conf = "./figures/fig.conf" 
    project_path = read_conf(conf,"project_path")
    results_sub_path = read_conf(conf,"results_sub_path");
    cor_data_file = read_conf(conf,"cor_data_file")
    annot_sub_path = read_conf(conf,"annot_sub_path")
    met_anno_fh = project_path*annot_sub_path*"ogm_meta_key_anno.csv"

    met_anno = CSV.read(met_anno_fh, DataFrame)
    cor_df = CSV.read(project_path*results_sub_path*cor_data_file,DataFrame)
    clean_met_anno!(cor_df,met_anno)


    # create_bargraphs
    sig = 0.05

    jhs_cor_sig_bar_col = gen_cor_sig_bar_col(cor_df,"jhs_p_qval","jhs_rho",sig)
    mesa_cor_sig_bar_col = gen_cor_sig_bar_col(cor_df,"mesa_p_qval","mesa_rho",sig)
    her_cor_sig_bar_col = gen_cor_sig_bar_col(cor_df,"her_p_qval","her_rho",sig)
    meta_cor_sig_bar_col = gen_cor_sig_bar_col(cor_df,"meta_p_qval","meta_rho",sig)

    stacked_bar_table = hcat(jhs_cor_sig_bar_col,mesa_cor_sig_bar_col,her_cor_sig_bar_col,meta_cor_sig_bar_col)
    stacked_bar_table = (stacked_bar_table / 1000)

    fig_cor_out = read_conf(conf,"fig_cor_out")
    bar_out_fh = project_path*fig_cor_out*"cor_bar.pdf"
    create_bargraphs(stacked_bar_table,bar_out_fh)


    meta_cor_sig_bar_col = gen_cor_sig_bar_col(cor_df,"meta_p_qval","meta_rho",sig)
    meta_bmi_cor_sig_bar_col = gen_cor_sig_bar_col(cor_df,"meta_p_bmi_qval","meta_rho_bmi",sig)
    meta_bmi_egfr_sig_bar_col = gen_cor_sig_bar_col(cor_df,"meta_p_bmi_egfr_qval","meta_rho_bmi_egfr",sig)

    adj_bar_table = vcat(sum(meta_cor_sig_bar_col[[1,3]]),sum(meta_bmi_cor_sig_bar_col[[1,3]]),sum(meta_bmi_egfr_sig_bar_col[[1,3]]))
    adj_bar_table = (adj_bar_table / 1000)
    bmi_pct = (adj_bar_table[2]/adj_bar_table[1] * 100)
    bmi_pct_str = @sprintf("%.0f",bmi_pct) * "%"
    egfr_pct = (adj_bar_table[3]/adj_bar_table[1] * 100)
    egfr_pct_str = @sprintf("%.0f",egfr_pct) * "%"

    adj_bar_out_fh = project_path*fig_cor_out*"adj_bar.pdf"
    create_adj_bargraphs(adj_bar_table, bmi_pct_str, egfr_pct_str, adj_bar_out_fh)

    # create_scatterplots

    jhs_mesa_scatter_df = gen_scatter_df(cor_df,"jhs_rho","mesa_rho")
    jhs_her_scatter_df = gen_scatter_df(cor_df,"jhs_rho","her_rho")
    mesa_her_scatter_df = gen_scatter_df(cor_df,"mesa_rho","her_rho")

    #create_scatter_plot(jhs_mesa_scatter_df,"JHS","MESA","jhs_mesa_rho_scatter.pdf",jhs_mesa_scatter_not_sig_df)
    jhs_mesa_scatter_fh = project_path*fig_cor_out*"jhs_mesa_rho_scatter.pdf"
    jhs_her_scatter_fh = project_path*fig_cor_out*"jhs_her_rho_scatter.pdf"
    mesa_her_scatter_fh = project_path*fig_cor_out*"mesa_her_rho_scatter.pdf"
    create_scatter_plot(jhs_mesa_scatter_df,"JHS Correlation Coefficient","MESA Correlation Coefficient",jhs_mesa_scatter_fh)
    create_scatter_plot(jhs_her_scatter_df,"JHS Correlation Coefficient","HERITAGE Correlation Coefficient",jhs_her_scatter_fh)
    create_scatter_plot(mesa_her_scatter_df,"MESA Correlation Coefficient","HERITAGE Correlation Coefficient",mesa_her_scatter_fh)

    # prep for create_heatmap
    met_order = CSV.read(project_path*annot_sub_path*"ordered_metabolites.csv",DataFrame)

    if size(met_order)[1] == size(met_anno)[1]
        println("USING SPECIFIED METABOLITE ORDER!!!")
        order_arr = reverse(sort_by(met_anno,"ogm_id",met_order))
        met_anno = copy(met_anno[order_arr,:])
    end

#    pm_class_order = ["Carbohydrates","Alkaloids","Nucleic acids","Combined_Lipids","Organoheterocyclic compounds","Benzenoids","Organic acids","Organic nitrogen compounds","Organic oxygen compounds"]
    #pm_class_order = ["Carbohydrates","Alkaloids","Nucleic acids","Glycerolipids","Glycerophospholipids","Sphingolipids","Fatty Acyls","Prenol Lipids","Sterol Lipids","Organoheterocyclic compounds","Benzenoids","Amino acids","Other Organic acids","Organic nitrogen compounds","Organic oxygen compounds"]
#    reverse!(pm_class_order)

#    met_anno[!,:pm_class_order] = copy(met_anno.pm_class)
#    for i in 1:length(pm_class_order)
#        t_rows = met_anno.pm_class .== pm_class_order[i]
#        met_anno.pm_class_order[t_rows] = "$(i)_".*met_anno.pm_class_order[t_rows]
#    end

    #==
    # split out and order lipid super classes
    lipid_order = ["Glycerolipids","Glycerophospholipids","Sphingolipids","Fatty Acyls","Prenol Lipids","Sterol Lipids"]
    for i in 1:length(lipid_order)
        lipid = lipid_order[i]

        for j in 1:size(met_anno)[1]
            if met_anno.RefMet_SuperClass[j] == lipid
                met_anno.pm_class_order[j] = met_anno.pm_class_order[j][1]*"$(i)"*met_anno.pm_class_order[j][2:end]
            end
        end

    end
    

    #split out AA from other organic acids
    aa_rows = met_anno.RefMet_SubClass .== "Amino acids"
    met_anno.pm_class_order[aa_rows] .= met_anno.pm_class_order[aa_rows][1][1] .* "1_Amino acids"

    ==#

    #met_anno[!,:RefMet_Main_Class_collapsed] = copy(met_anno.RefMet_Main_Class)
    #for main in met_anno.RefMet_Main_Class
    #    #t_df = filter(row -> row.RefMet_Main_Class == main, met_anno)
    #    t_rows = met_anno.RefMet_Main_Class .== main
    #    if sum(t_rows) < 10
    #        met_anno.RefMet_Main_Class_collapsed[t_rows] = met_anno.pm_class[t_rows]
    #    end
    #end

#    sort!(met_anno,[:pm_class_order,:pm_subclass])
    met_anno = redo_bw_colors(met_anno)

    (cor_mat, znorm_cor_mat, cor_mat_sig, p_mat, q_mat, prots, mets) = gen_cor_mat(cor_df,met_anno)
    t_cor_mat = trim_cor_mat(cor_mat,0.3)
    t_znorm_cor_mat = trim_cor_mat(znorm_cor_mat,3)

    # reorder cor_mat
    t_fh = project_path*fig_cor_out*"t_heatmap.pdf"
    t_hmap = create_heatmap(t_cor_mat,met_anno,t_fh,"NA","NULL")
    prot_order = t_hmap[Symbol("colInd")]
    met_order = t_hmap[Symbol("rowInd")]
    t_cor_mat2 = copy(t_cor_mat[prot_order,:])

    #==
    pm_class_order2 = sort(collect(Set(met_anno.pm_class_order)))
    met_order = []
    jump = 0
    #for i in 1:length(pm_class_order)
    for i in 1:length(pm_class_order2)
        println(i)
        #class = pm_class_order[i]
        class = pm_class_order2[i]
        main = ""

        #for main in sort(collect(Set(filter(row -> row.pm_class == class, met_anno).RefMet_Main_Class_collapsed)))
            #c_rows = (met_anno.pm_class .== class)# .* (met_anno.RefMet_Main_Class_collapsed .== main)
            c_rows = (met_anno.pm_class_order .== class)
            c_t_cor_mat = t_cor_mat2[:,c_rows]
            c_met_anno = met_anno[c_rows,:]

            if size(c_t_cor_mat)[2] > 1
                t_fh2 = project_path*fig_cor_out*"t_heatmap2_$(class)_$(main).pdf"
                t_hmap2 = create_heatmap(c_t_cor_mat,c_met_anno,t_fh2,"NULL","NA")

                met_order = vcat(met_order,t_hmap2[Symbol("rowInd")].+jump)
            else
                met_order = vcat(met_order,[1+jump])
            end
            jump = jump + sum(c_rows)
        #end
    end
    met_order = Int.(met_order)
    t_cor_mat3 = copy(t_cor_mat2[:,met_order])
    met_anno3 = copy(met_anno[met_order,:])

    t_fh3 = project_path*fig_cor_out*"t_heatmap3.pdf"
    t_hmap3 = create_heatmap(t_cor_mat3,met_anno3,t_fh3,"NA","NA")


    # reorder cor_mat z_norm
    t_fh = project_path*fig_cor_out*"t_heatmap_znorm.pdf"
    t_hmap = create_heatmap(t_cor_mat,met_anno,t_fh,"NULL","NULL")
    prot_order_znorm = t_hmap[Symbol("colInd")]
    t_cor_mat2 = copy(t_znorm_cor_mat[prot_order_znorm,:])

    met_order = []
    jump = 0
    for i in 1:length(pm_class_order)
        println(i)
        class = pm_class_order[i]

        c_rows = met_anno.pm_class .== class
        c_t_cor_mat = t_cor_mat2[:,c_rows]
        c_met_anno = met_anno[c_rows,:]

        t_fh2 = project_path*fig_cor_out*"t_heatmap2_znorm_$(class).pdf"
        t_hmap2 = create_heatmap(c_t_cor_mat,c_met_anno,t_fh2,"NULL","NA")

        met_order = vcat(met_order,t_hmap2[Symbol("rowInd")].+jump)
        jump = jump + sum(c_rows)
    end
    met_order_znorm = Int.(met_order)
    t_cor_mat3 = copy(t_cor_mat2[:,met_order])
    met_anno3_znorm = copy(met_anno[met_order_znorm,:])

    t_fh3 = project_path*fig_cor_out*"t_heatmap3_znorm.pdf"
    t_hmap3 = create_heatmap(t_cor_mat3,met_anno3_znorm,t_fh3,"NA","NA")

    ##==#

    t_fh3 = project_path*fig_cor_out*"t_heatmap3.pdf"
    t_hmap2 = create_heatmap(t_cor_mat,met_anno,t_fh3,"NULL",t_hmap[Symbol("colInd")])

    # create_heatmap
    fh = project_path*fig_cor_out*"heatmap_znorm_0.pdf"
    hmap = create_heatmap(t_cor_mat,met_anno,fh,"NA","NULL")

    fh_znorm = project_path*fig_cor_out*"heatmap_znorm_1.pdf"
    hmap_znorm = create_heatmap(t_znorm_cor_mat,met_anno,fh_znorm,"NA","NULL")

    #==#

    prot_anno = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_source_data/annotations/kegg_annotated_proteins.csv", DataFrame)
    sort!(prot_anno,:soma_id)
    filter!(row -> in(row.soma_id,Set(prots)), prot_anno)

    clustered_prot_df = prot_anno[hmap[:colInd],[:soma_id,:target_full_name,:target]]
    CSV.write(project_path*results_sub_path*"hmap_prots_order.csv",clustered_prot_df)

    apoe_sorted_mets = sort(filter(row -> row.PROT == "SL000276",cor_df),:meta_p_qval).OGM
    for met in apoe_sorted_mets
        t_met_anno = filter(row -> row.ogm_id == met, met_anno)
        println(t_met_anno.pm_class[1])
    end
    ==#

    ### hormone correlation examples
    # Insulin	SL000021
    # Ghrelin	SL004269
    # Adiponectin	SL004258
    # FGF19	SL004337

    hormones = ["SL000021","SL004269","SL004258","SL004337"]

    for hormone in hormones
        #hormone = hormones[1]
        println(hormone)
        h = prots .== hormone
        h_arr = sign.(cor_mat[h,:]) .* -log.(10,p_mat[h,:])
        h_cor_arr = cor_mat[h,:]
        #h_arr = cor_mat[h,:]
        #h_arr = reverse!(copy(h_arr[met_order]))
        h_arr = copy(h_arr[1:length(h_arr)])
        #h_arr = reverse!(copy(h_arr[1:length(h_arr)]))
        #bar(h_arr)

        #bar_colors = []
        #for h in h_arr
        #    if h > 0
        #        push!(bar_colors,"firebrick3")
        #    else
        #        push!(bar_colors,"navy")
        #    end
        #end
        #bar_colors_df = DataFrame(colors = bar_colors)

        out_fh = project_path*fig_cor_out*"hormone_bar_$(hormone).pdf"
        trim = 0.3
        @rput(h_arr,h_cor_arr,out_fh,trim)
        R"""
            pdf(out_fh, width = 3, height = 10)
            par(mar=c(0, 0, 0, 0))

            bar_col = sign(h_arr)
            bar_col[bar_col == -1] = "navy"
            bar_col[bar_col == 1] = "firebrick3"

            cols = colorRampPalette(c("navy", "white", "firebrick3"))(15)
            ranges = seq(-trim,trim,length.out = length(cols)-1)
            ranges = c(-1,ranges,1)

            bar_col = c()
            for (i in 1:length(h_arr))
            {
                for (j in 1:(length(ranges)-1))
                {
                    if (h_cor_arr[i] > ranges[j] && h_cor_arr[i] <= ranges[j+1])
                    {
                        bar_col = c(bar_col,cols[j])
                    }
                }
            }


            ymax = 1.2*max(abs(h_arr))
            barplot(h_arr, col = bar_col, border = NA, xlim = c(-ymax,ymax), horiz = TRUE, axes = FALSE)

            #barplot(h_arr, col = c("navy","firebrick3"), ylab = expression(paste("Correlation Tests (n x ",10^{3},")")), ylim = c(0,500), cex.names = 2, cex.lab = 2, cex.axis = 2, legend.text = c("Negative Correlation", "No Correlation (q > 0.05)", "Positive Correlation"), args.legend = list(x = "bottom", inset = c(0,-.35), cex = 2))

            dev.off()
        """
    end

end

do_fig_cor()









