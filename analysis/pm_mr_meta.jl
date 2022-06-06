using DataFrames,CSV,RCall

mr = CSV.read("/Users/aaroneisman/Google Drive/AMS/BCBI/projects/benson/OG_papers/protein-metabolite/pm_paper/pm_data/pm_results_data/mr/20220429_all_pm_mr.csv", DataFrame, missingstring = "null")

mr_methods = ["ivw","maxlik","med","med_wt","egg","egg_int"]
#mr_methods = ["maxlik","med","med_wt","egg","egg_int"]

function call_metafor(beta,se)
    @rput(beta,se)
    R"""
        options(warn=0)
        require('metafor')
        m = rma(beta, se^2, method="FE",verbose=FALSE);
    """
    @rget(m)

    return m
end

for mr_method in mr_methods
    println("Starting $(mr_method) meta-analysis...")
    #mr_method = "ivw"

    mr[!,Symbol("META_B_"*mr_method)] .= 999.99
    mr[!,Symbol("META_SE_"*mr_method)] .= 999.99
    mr[!,Symbol("META_P_"*mr_method)] .= 999.99

    for i in 1:size(mr)[1]
        beta_i = [mr[i,Symbol("r.JHS_B_"*mr_method)],mr[i,Symbol("r.MESA_B_"*mr_method)],mr[i,Symbol("r.HERITAGE_B_"*mr_method)]]
        se_i = [mr[i,Symbol("r.JHS_SE_"*mr_method)],mr[i,Symbol("r.MESA_SE_"*mr_method)],mr[i,Symbol("r.HERITAGE_SE_"*mr_method)]]

        if sum(ismissing.(beta_i)) < 3 && sum(ismissing.(se_i)) < 3
            m_i = call_metafor(beta_i,se_i)

            mr[i,Symbol("META_B_"*mr_method)] = m_i[:beta]
            mr[i,Symbol("META_SE_"*mr_method)] = m_i[:se]
            mr[i,Symbol("META_P_"*mr_method)] = m_i[:pval]
        end

        if (i % 10000) == 0
            println(i)
        end

    end

    println("...Finished!")

end