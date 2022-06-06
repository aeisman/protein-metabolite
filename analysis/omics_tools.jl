using DataFrames, Statistics, GLM, HypothesisTests, Printf, LinearAlgebra, Clustering

function read_conf(fh,var)
    
    var_val = ""
    for line in readlines(fh)
        line_arr = split(line,",")
        if line_arr[1] == var
            var_val = line_arr[2]
        end
    end

    if var_val == ""
        @error("$(var) not in $(fh)")
    end
    
    return var_val
end

###################################################
### TOOLS TO NORMALIZE AND TRANSFORM OMICS DATA ###
###################################################

function ln_scale(arr,mu,sd)
    ### Function ln transforms and scales an array of data to chosen mu and sd.
    # returns a transformed array of the same length
    # ignores missing values for calculation but retains them in returned array to preserve index numbers
    t_arr = arr
    t_arr_ln = log.(t_arr)
    sd_raw = Statistics.std(collect(skipmissing(t_arr_ln))) ## skip missing values for calculating standard deviation
    mu_raw = Statistics.mean(collect(skipmissing(t_arr_ln))) ## skip missing values for calculating mean
    
    t_arr_ln_scale = (t_arr_ln .+ (mu-mu_raw)) .* sd/sd_raw ## add difference between desired mu and mu_raw and multiply by ratio of desired to raw sd
    
    return t_arr_ln_scale
end

function by_batch!(var,df,mu,sd,b_var)
    ### Function to call ln_scale by batch according to a batch variable.
    # Tranformations are applied to the provided dataframe
    t_df = dropmissing(df[:,[b_var,var]])
    #batches = Set(df[:,b_var])
    batches = Set(t_df[:,b_var])
    #batches = setdiff(batches,Set([missing]))
    for batch in batches
        #println("Batch: $batch")
        #println("ID batch")
        in_batch = df[:,b_var] .== batch
        #println("Edit in_batch")
        if length(in_batch[ismissing.(in_batch)]) > 0
            in_batch[ismissing.(in_batch)] .= false
            in_batch = Bool.(in_batch)
        end
        #println("Create batch arr")
        batch_arr = df[in_batch,var]
        #println("Generate transformed array")
        t_batch_arr = ln_scale(batch_arr,mu,sd)
        #println("Assign transformed array")
        df[in_batch,var] = t_batch_arr
        #println("Transformed array")
    end
    return df
end

function age_sex_adj!(var,df,age_var,sex_var)
    ### Function to replace a variable in a dataframe by the residuals of var ~ age + sex + res
    not_missing = .!ismissing.(df[:,var])
    df_temp = DataFrame()
    df_temp.y = df[:,var]
    df_temp.age = df[:,age_var]
    df_temp.sex = categorical(df[:,sex_var])
    mod = lm(@formula(y ~ age + sex), df_temp)
    res = residuals(mod)
    df[not_missing,var] = res
    return df
end

function var_adj!(var,df,adj_var)
    ### Function to replace a variable in a dataframe by the residuals of var ~ adj_var + res
    not_missing = .!ismissing.(df[:,var])
    adj_missing = ismissing.(df[:,adj_var])
    make_missing = not_missing .* adj_missing
    complete = not_missing .* .!adj_missing
    if sum(make_missing) > 0
        df[make_missing,var] .= missing
    end
    df_temp = DataFrame()
    df_temp.y = df[:,var]
    df_temp.x = df[:,adj_var]
    mod = lm(@formula(y ~ x), df_temp)
    res = residuals(mod)
    df[complete,var] = res
    return df
end

function var_adj_2!(var_arr,adj_var_arr)
    ### Function to take an array of numbers and replace with residuals of var ~ adj_var + res
    not_missing = .!ismissing.(df[:,var])
    df_temp = DataFrame()
    df_temp.y = df[:,var]
    df_temp.x = df[:,adj_var]
    mod = lm(@formula(y ~ x), df_temp)
    res = residuals(mod)
    df[not_missing,var] = res
    return df
end

function adj_vars!(var_set,df,mu,sd,b_var,adj_var)
    ### Wrapper function to adjust a set of vars (ln scale and then age/sex adjust) by batch
    #b_var == 999 ==> do not log transform and batch adjust (already done)
    for var in var_set
        println(var)
        var_symb = Symbol(var)
        if b_var != 999 ### treat as single batch
            df = by_batch!(var_symb,df,mu,sd,b_var)
        end
        df = var_adj!(var_symb,df,Symbol(adj_var))
    end
    return df
end

################################
### OMICS ANALYSIS FUNCTIONS ###
################################

function pearson_cor_test(a,b)
    ### Function to calculate the pearson correlation between two arrays.
    # Return pearson correlation and associated Fisher transformed p-value and number of paird values
    # Returns a correlation of 0 with a p-value of 999 if length of arrays after dropping missing values is less than 5
    # Drops incomplete pairs (e.g. if a value is missing for a, b, or both for a given row)

    if length(a) != length(b)
        error("Arrays provided are different lengths")
    end
    df_temp = DataFrame()
    df_temp.a = a
    df_temp.b = b
    df_temp = dropmissing(df_temp)
    if size(df_temp)[1] < 10
        return([0,999,0])
    end
    rho = cor(df_temp.a,df_temp.b)
    ### Fisher Transformation to calculated p-value ###
    ### gives slightly different values than R at small length (likely T vs Z test / degrees of freedom)
    p = pvalue(OneSampleZTest(atanh(rho),1,length(df_temp.a)-3))
    return(rho,p,length(df_temp.a))
end

function add_classes!(cor_df, var_class_key)
    cor_df.var1_class = ""
    cor_df.var2_class = ""

    for i in 1:size(cor_df)[1]
        if haskey(var_class_key, cor_df.var1[i])
            cor_df.var1_class[i] = var_class_key[cor_df.var1[i]]
        end
        if haskey(var_class_key, cor_df.var2[i])
            cor_df.var2_class[i] = var_class_key[cor_df.var2[i]]
        end
    end

    return cor_df
end

function add_pqtl!(cor_df, pqtl_set)
    jhs_prot_met_cor_df[:,:var1_pqtl] .= 0
    jhs_prot_met_cor_df[:,:var2_pqtl] .= 0

    for i in 1:size(cor_df)[1]
        if in(cor_df.var1[i],pqtl_set)
            cor_df.var1_pqtl[i] = 1
        end
        if in(cor_df.var2[i],pqtl_set)
            cor_df.var2_pqtl[i] = 1
        end
    end

    return cor_df
end

function find_correlations(vars1,vars2,df,var_key)
    ### Wrapper function to find pearson correlations for all combinations of vars1 and vars2.
    # Returns correlations as a dataframe
    # Will also add alternate names for the variables if available within the Dict var_key

    cor_results = DataFrame(var1 = String[], var1_alt = String[], var2 = String[], var2_alt = String[], rho = Float64[], p = Float64[], n = Int[])
    total = length(vars1)*length(vars2)
    i = 0
    for var1 in vars1
        for var2 in vars2
            #println(var2)
            if var1 != var2
                #println("$var1 $var2")
                if rem(i,round(0.05*total)) == 0
                    pct = i/total*100
                    pct_str = @sprintf("%.1f",pct)
                    println("$i of $total | $pct_str%")
                end
                var1_alt = ""
                var2_alt = ""
                if haskey(var_key,String(var1))
                    var1_alt = var_key[String(var1)]
                end
                if haskey(var_key,String(var2))
                    var2_alt = var_key[String(var2)]
                end
                (rho,p,n) = pearson_cor_test(df[:,Symbol(var1)],df[:,Symbol(var2)])
                push!(cor_results,[String(var1),var1_alt,String(var2),var2_alt,rho,p,n])
                i += 1
            end
        end
    end
    return(cor_results)
end

function add_adj_correlations!(vars1,vars2,df,cor_df,var_key,adj_var,adj_var2)
    temp_df = copy(df,copycols = true)
    vars_to_adj = collect(union(Set(vars1),Set(vars2)))
    cor_name_sym = Symbol("rho_"*string(adj_var))
    p_name_sym = Symbol("p_"*string(adj_var))

    adj_vars!(vars_to_adj,temp_df,0,1,999,Symbol(adj_var))
    if adj_var2 != 999 ### only adjust with one variable
        adj_vars!(vars_to_adj,temp_df,0,1,999,Symbol(adj_var2))
        cor_name_sym = Symbol("rho_"*string(adj_var)*"_"*string(adj_var2))
        p_name_sym = Symbol("p_"*string(adj_var)*"_"*string(adj_var2))
    end
    temp_cor_df = find_correlations(vars1,vars2,temp_df,var_key)

    ## check that cor_df are in the same order
        if cor_df.var1 != temp_cor_df.var1 || cor_df.var2 != temp_cor_df.var2
            error("Correlation dataframe variables are in different orders. Results will not line up!")
        end
    ##

    cor_df[:,cor_name_sym] = temp_cor_df.rho
    cor_df[:,p_name_sym] = temp_cor_df.p

    return cor_df
end


### COLOR PALETTE FUNCTIONS ###
#get colors using R library // COME BACK TO ONLY INCLUDE ONES WE NEED
function get_viridis_colors(n)
    @rput(n)
    R"""
    library(viridis)
        colors = viridis(n)
    """
    @rget(colors)
    return colors
end

function get_unikn_colors(n,pal)

    @rput(n,pal)
    R"""
        library(unikn)

        if (pal == 1) {
            colors = usecol(pal_seegruen,n)
        } else if (pal == 2) {
            colors = usecol(pal_bordeaux,n)
        } else if (pal == 3) {
            colors = usecol(pal_seeblau,n)
        } else {
            colors = usecol(pal_seegreen,n)
        }

    """
    @rget(colors)

    return colors
end

function get_brewer_pal_colors(n,name)

    @rput(n,name)
    R"""
        library(RColorBrewer)
        colors = brewer.pal(n,name)
    """
    @rget(colors)

    return colors
end

function get_colorRampPalette(col_arr, n)
    @rput(col_arr,n)
    R"""
        colors = colorRampPalette(col_arr)(n)
    """
    @rget(colors)

    return colors
end