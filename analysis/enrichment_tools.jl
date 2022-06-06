using CSV, DataFrames, RCall, Distributions, Plots, StatsBase, Random

### EA Function
# take data table, var1, var2, var1_cat, var2_cat, use_cor_wt (0 or 1)
function ogea(arr,p)
    ### test parameters
    #p = 1
    # add rank column to arr
    arr = hcat(arr,1:size(arr)[1])

    S = arr[arr[:,1] .== 1,:]
    not_S = arr[arr[:,1] .== 0,:]

    n = size(arr)[1]
    n_hits = size(S)[1]


    N_R = 0
    for j in 1:size(S)[1]
        N_R += abs(S[j,2])^p
    end

    cum_sum = []
    for i in 1:n
        P_hit = 0
        P_miss = 0
        for j in 1:n_hits
            if S[j,3] <= i
                P_hit += abs(S[j,2])^p / N_R
                #P_hit += abs(S[j,2])^p / N_R * (n - n_hits)
            end
        end

        for j in 1:size(not_S)[1]
            if not_S[j,3] <= i
                P_miss += 1 / (n - n_hits)
                #P_miss += 1
            end
        end

        push!(cum_sum,P_hit-P_miss)
    end

    return cum_sum
end


function ogea_es(ogea_cum_sum)
    ## Decided only positive enrichments make sense for OGEA -- this is 2/2 the choice to sort by p-value instead of by correlation coefficient
    es = maximum(ogea_cum_sum)
    #min_es = minimum(ogea_cum_sum)
    #if abs(min_es) > es
    #    es = min_es
    #end

    return es
end

function ogea_sig_thresh(n_pos,n_neg,n_itr,nom,cor_data,p,assume_normal,sortby,var2_distr)
    # function to simulate and determine significance threshold for a given n_pos, n_neg using random random walk

    n = n_pos + n_neg

    if p == 0
        println("Using p = 0")
        es_arr = []
        for i in 1:n_itr
            rand_list = rand(n)
            rand_list_bin = rand_list .< n_pos / n
            rand_list_bin_cor = hcat(rand_list_bin,ones(length(rand_list_bin)))
            ogea_cum_sum = ogea(rand_list_bin_cor,p)
            es = ogea_es(ogea_cum_sum)
            push!(es_arr,es)
        end
        sig_h = sort(es_arr)[convert(Int,ceil((1-nom)*n_itr))]
        sig_l = sort(es_arr)[convert(Int,floor((nom)*n_itr))]

    else
        println("Using p = $(p)")
        all_cor = cor_data.var_cor
        mu = mean(all_cor)
        sd = std(all_cor)

        pos_cor = filter(row -> row.var2_cat == true, cor_data).var_cor
        mu_pos = mean(pos_cor)
        sd_pos = std(pos_cor)
        
        neg_cor = filter(row -> row.var2_cat == false, cor_data).var_cor
        mu_neg = mean(neg_cor)
        sd_neg = std(neg_cor)

        cor_dist = Normal(mu,sd) ## assume normal distribution of rho w/ mean mu and sigma of sd
        cor_pos_dist = Normal(mu_pos,sd_pos)
        cor_neg_dist = Normal(mu_neg,sd_neg)

        es_arr = []
        for i in 1:n_itr
            #rand_list = rand(n)
            #rand_list_bin = rand_list .< n_pos / n

            cat_list = vcat(ones(n_pos),zeros(n_neg))
            
            if assume_normal == 1
                if var2_distr == 1
                    #new way
                    rand_cor_pos = rand(cor_pos_dist,n_pos)
                    rand_cor_neg = rand(cor_neg_dist,n_neg)
                    rand_cor = vcat(rand_cor_pos,rand_cor_neg)
                else
                    rand_cor = rand(cor_dist,n)
                end
            else
                if var2_distr == 1
                    rand_cor_pos = shuffle(pos_cor)[1:n_pos]
                    rand_cor_neg = shuffle(neg_cor)[1:n_neg]
                    rand_cor = vcat(rand_cor_pos,rand_cor_neg)
                else
                    rand_cor = shuffle(all_cor)[1:n]
                end
            end

            #rand_list_bin_cor = hcat(rand_list_bin,rand_cor)

            # new way
            if sortby == "abs"
                rand_cor_perm = sortperm(rand_cor, by = abs, rev = true)
            else
                rand_cor_perm = sortperm(rand_cor, rev = true)
            end
            rand_list_bin_cor = hcat(cat_list[rand_cor_perm],rand_cor[rand_cor_perm])

            ogea_cum_sum = ogea(rand_list_bin_cor,p)
            es = ogea_es(ogea_cum_sum)
            push!(es_arr,es)
        end
        es_arr_pos = es_arr[es_arr .> 0]
        es_arr_neg = es_arr[es_arr .< 0]

        #don't separate by positive and negative es if using absolute value sorting
        if sortby == "abs"
            es_arr_pos = es_arr
            es_arr_neg = es_arr
        end
        
        if convert(Int,ceil((1-nom)*length(es_arr_pos))) == 0
            sig_h = 999999999
        else
            sig_h = sort(es_arr_pos)[convert(Int,ceil((1-nom)*length(es_arr_pos)))]
        end
        if convert(Int,floor((nom)*length(es_arr_neg))) == 0
            sig_l = -999999999
        else
            sig_l = sort(es_arr_neg)[convert(Int,floor((nom)*length(es_arr_neg)))]
        end
    end

    return (sig_h,sig_l,es_arr)
end

function ogea_graph(ogea_graph_df,vars,fn,labels,colors,incl_legend)
    #vars = ["SL003650","SL005574"]
    # generate array of maximums to plot points
    maxes = []
    idxs = []
    for var in vars
        max_temp = findmax(ogea_graph_df[:,Symbol(var)])
        min_temp = findmin(ogea_graph_df[:,Symbol(var)])
        #if abs(min_temp[1]) > abs(max_temp[1])
        #    max_temp = min_temp
        #end
        push!(maxes,max_temp[1])
        push!(idxs,max_temp[2])
    end
    max_idx = hcat(idxs,maxes)
    m = maximum(maxes)


    @rput(ogea_graph_df,vars,max_idx,m,fn,labels,colors,incl_legend)
    R"""
        library(dplyr)
        library(viridis)
        library(unikn)

        legend_n = min(c(length(labels),10))
        
        
        if (incl_legend == 0){
            pdf(fn, width = 14, height = 7)
            par(mar=c(2, 3, 1, 1),pin = c(6.75,2.5))
        } else {
            pdf(fn, width = 9.5, height = 7)
            par(mar=c(5, 5, 1, 12),pin = c(4.75,4.5))
        }
        

        n_mets = dim(ogea_graph_df)[1]

        #colors = usecol(pal_seegruen,n = length(vars))
        #colors = viridis(length(vars))

        matplot(select(ogea_graph_df,as.character(vars)),type = "l", ylab = "Enrichment Score (ES)", xlab = "Ranked Metabolites", xaxs = "i", xaxt='n', yaxt='n', xlim = c(0,n_mets), col = colors, lwd = 2, cex.lab = 2)
        points(max_idx,col = colors, pch = 16, cex = 2)
        abline(h=0)
        par(xpd = TRUE)

        axis(tick = TRUE, side = 1, at = c(0,n_mets), cex.axis = 1.5)

        #text(n_mets,min(ogea_graph_df),labels = c(as.character(n_mets)))
        #mtext(n_mets, side = 1, adj = 1, cex = 1.5, line = 0.5)
        
        if (incl_legend == 1) {
            legend("topleft", inset = c(1.03,0),   # Coordinates (x also accepts keywords)
                as.character(labels[1:legend_n]), # Vector with the name of each group
                fill = colors, #viridis(length(vars)),   # Creates boxes in the legend with the specified colors
                #col = viridis(length(vars)), # Color of lines or symbols
                border = "black", # Fill box border color
                ncol = 1,
                #lty = 1, lwd = 1,         # Line type and width
                #pch = 1,              # Add pch symbols to legend lines or boxes
                #bty = "o",        # Box type (bty = "n" removes the box)
                #bg = par("bg")    # Background color of the legend
                #box.lwd = par("lwd"), # Legend box line width
                #box.lty = par("lty"), # Legend box line type
                #box.col = par("fg"),  # Legend box line color
                cex = 1.5,          # Legend size
                #horiz = FALSE     # Horizontal (TRUE) or vertical (FALSE) legend
                #title = NULL      # Legend title
            )
        }

        dev.off() 
    """
end

function ogea_es_graph(ogea_results,fn)
    
    sig_n = sum(ogea_results.sig .== true)
    
    @rput(ogea_results,fn,sig_n)
    R"""
    library(dplyr)
    library(viridis)
    
    pdf(fn,width = 10, height = 3)

    par(mar=c(5, 5, 1, 1.5))

    n_prots = length(ogea_results$es)

    plot(as.numeric(ogea_results$es), type = "l", ylab = "ES", xlab = "Ranked Proteins", xaxs = "i", xaxt='n', lwd = 2, yaxt = 'n', cex.lab = 2, xlim = c(0,n_prots))
    abline(v=sig_n,col="red",lwd = 2)

    axis(tick = TRUE, side = 1, at = c(0,sig_n,n_prots), cex.axis = 2)

    dev.off() 
"""
end

function calc_p_from_dist(dist,val)

    n_greater = 0
    for i in 1:length(dist)
        if dist[i] > val
            n_greater += 1
        end
    end

    p_est = n_greater / length(dist)
    if p_est == 0
        p_est = 1 / length(dist)
    end

    return p_est
end