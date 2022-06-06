using DataFrames,CSV
include("../analysis/omics_tools.jl")

function create_ogm_meta(ogm_key_fh, ogm_studies_arr, ogm_platforms_arr)
    ## Platform priority: C8 --> HILIC --> AMIDE

    ogm_key = CSV.read(ogm_key_fh,DataFrame)
    filter!(row -> !ismissing(row.ogm_name),ogm_key)
    ogm_name_set = Set(ogm_key.ogm_name)
    
    studies = ogm_studies_arr

    ogm_meta_key = DataFrame()

    for ogm_name_temp in ogm_name_set
        println(ogm_name_temp)

        ogm_name_subdf = filter(row -> row.ogm_name == ogm_name_temp, ogm_key)

        platform_set = intersect(Set(ogm_name_subdf.ogm_platform),Set(ogm_platforms_arr))

        if in("c8",platform_set)
            println("Choose c8!")

            row = filter(row -> row.ogm_platform == "c8", ogm_name_subdf)

            for study in studies
                study_name = Symbol(study*"_name")
                study_id = Symbol(study*"_id")
                study_platform = Symbol(study*"_platform")

                study_name_c8 = !ismissing(row[:,Symbol(study_name)])
                study_id_c8 = !ismissing(row[:,Symbol(study_id)])

                (study_name_hilic,study_id_hilic) = (false,false)
                if in("hilic",platform_set)
                    study_name_hilic = !ismissing(filter(row -> row.ogm_platform == "hilic", ogm_name_subdf)[:,study_name][1])
                    study_id_hilic = !ismissing(filter(row -> row.ogm_platform == "hilic", ogm_name_subdf)[:,study_id][1])
                end

                (study_name_amide,study_id_amide) = (false,false)
                if in("amide",platform_set)
                    study_name_amide = !ismissing(filter(row -> row.ogm_platform == "amide", ogm_name_subdf)[:,study_name][1])
                    study_id_amide = !ismissing(filter(row -> row.ogm_platform == "amide", ogm_name_subdf)[:,study_id][1])
                end

                if study_name_c8 || study_id_c8
                    #keep row as is
                elseif in("hilic",platform_set) && (study_name_hilic || study_id_hilic)
                    row[:,study_name] = ogm_name_subdf[ogm_name_subdf.ogm_platform .== "hilic", study_name]
                    row[:,study_id] = ogm_name_subdf[ogm_name_subdf.ogm_platform .== "hilic", study_id]
                    row[:,study_platform] = ogm_name_subdf[ogm_name_subdf.ogm_platform .== "hilic", study_platform]
                elseif in("amide",platform_set) && (study_name_amide || study_id_amide)
                    row[:,study_name] = ogm_name_subdf[ogm_name_subdf.ogm_platform .== "amide", study_name]
                    row[:,study_id] = ogm_name_subdf[ogm_name_subdf.ogm_platform .== "amide", study_id]
                    row[:,study_platform] = ogm_name_subdf[ogm_name_subdf.ogm_platform .== "amide", study_platform]
                end
            end
            ogm_meta_key = vcat(ogm_meta_key,row)

        elseif in("hilic",platform_set)
            println("Choose hilic!")

            row = filter(row -> row.ogm_platform == "hilic", ogm_name_subdf)

            for study in studies
                study_name = Symbol(study*"_name")
                study_id = Symbol(study*"_id")
                study_platform = Symbol(study*"_platform")

                study_name_hilic = !ismissing(filter(row -> row.ogm_platform == "hilic", ogm_name_subdf)[:,study_name][1])
                study_id_hilic = !ismissing(filter(row -> row.ogm_platform == "hilic", ogm_name_subdf)[:,study_id][1])

                (study_name_amide,study_id_amide) = (false,false)
                if in("amide",platform_set)
                    study_name_amide = !ismissing(filter(row -> row.ogm_platform == "amide", ogm_name_subdf)[:,study_name][1])
                    study_id_amide = !ismissing(filter(row -> row.ogm_platform == "amide", ogm_name_subdf)[:,study_id][1])
                end

                if study_name_hilic || study_id_hilic
                    #keep row as is
                elseif (in("amide",platform_set)) && study_name_amide || study_id_amide
                    row[:,study_name] = ogm_name_subdf[ogm_name_subdf.ogm_platform .== "amide", study_name]
                    row[:,study_id] = ogm_name_subdf[ogm_name_subdf.ogm_platform .== "amide", study_id]
                    row[:,study_platform] = ogm_name_subdf[ogm_name_subdf.ogm_platform .== "amide", study_platform]
                end
            end

            ogm_meta_key = vcat(ogm_meta_key,row)

        elseif in("amide",platform_set)
            println("Choose amide!")

            row = filter(row -> row.ogm_platform == "amide", ogm_name_subdf)
            
            ogm_meta_key = vcat(ogm_meta_key,row)
        end
        println()
    end

    return(ogm_meta_key)

end

function gen_met_anno(met_anno_fh)
    met_anno = CSV.read(met_anno_fh,DataFrame)
    filter!(row -> !ismissing(row.RefMet_SuperClass),met_anno)
    sort!(met_anno,[:RefMet_SuperClass,:RefMet_Main_Class])

    pm_lipids = Set(["Fatty Acyls","Glycerolipids","Glycerophospholipids","Sphingolipids","Sterol Lipids","Prenol Lipids","Polyketides","Saccharolipids"])


    met_anno[!,:pm_class] .= ""
    for i in 1:size(met_anno)[1]
        if in(met_anno[i,:RefMet_SuperClass],pm_lipids)
            met_anno[i,:pm_class] = "Combined_Lipids"
        elseif ismissing(met_anno[i,:RefMet_SuperClass])
            met_anno[i,:pm_class] = ""
            met_anno[i,:RefMet_SuperClass] = ""
            met_anno[i,:RefMet_Main_Class] = ""
            met_anno[i,:RefMet_SubClass] = ""
        else
            met_anno[i,:pm_class] = met_anno[i,:RefMet_SuperClass]
        end
    end

    met_anno[!,:pm_subclass] .= ""
    for i in 1:size(met_anno)[1]
        if met_anno.pm_class[i] == "Combined_Lipids"
            met_anno.pm_subclass[i] = met_anno.RefMet_SuperClass[i]
        elseif met_anno.RefMet_SubClass[i] == "Amino acids"
            met_anno.pm_subclass[i] = met_anno.RefMet_SubClass[i]
        else
            met_anno.pm_subclass[i] = met_anno.RefMet_Main_Class[i]
        end
    end

    sort!(met_anno,[:pm_class,:RefMet_SuperClass,:RefMet_Main_Class,:RefMet_SubClass])

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

conf = "./data_prep/data_prep.conf"
project_path = read_conf(conf,"project_path")
ogm_sub_path = read_conf(conf,"ogm_sub_path")
ogm_key_fn = read_conf(conf,"ogm_key_fn")
ogm_key_fh = project_path*ogm_sub_path*ogm_key_fn
ogm_studies_arr = split(read_conf(conf,"ogm_studies"),"|")
ogm_platforms_arr = split(read_conf(conf,"ogm_platforms"),"|")

ogm_meta_key = create_ogm_meta(ogm_key_fh, ogm_studies_arr, ogm_platforms_arr)

CSV.write(project_path*ogm_sub_path*"ogm_meta_key.csv",ogm_meta_key)

met_anno = gen_met_anno(project_path*ogm_sub_path*"ogm_meta_key.csv")

CSV.write(project_path*ogm_sub_path*"ogm_meta_key_anno.csv",met_anno)