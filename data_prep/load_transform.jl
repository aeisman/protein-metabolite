
function load_transform(studies_arr)
    for study in studies_arr
        load_transform_script = "../data_prep/"*study*"_load_transform.jl"
        include(load_transform_script)
    end
end

studies_arr = ["jhs","mesa","heritage"]
load_transform(studies_arr)