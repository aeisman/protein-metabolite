using DataFrames,CSV,StatsBase,RCall
try
    include("../analysis/omics_tools.jl")
catch
    include("./analysis/omics_tools.jl")
end

### Load data paths from data_prep.conf config file
conf = "./data_prep/data_prep.conf"
project_path = read_conf(conf,"project_path")
mesa_data_sub_path = read_conf(conf,"mesa_data_sub_path")

mesa_clin_fn = read_conf(conf,"mesa_clin_fn")
mesa_linker_fn = read_conf(conf,"mesa_linker_fn")
mesa_hil_fn = read_conf(conf,"mesa_hil_fn")
mesa_ami_fn = read_conf(conf,"mesa_ami_fn")
mesa_prot_fn = read_conf(conf,"mesa_prot_fn")

### Join all MESA files into single table

### load clinical
clin_fh = project_path*mesa_data_sub_path*mesa_clin_fn
clin_df = CSV.read(clin_fh,DataFrame)
for i in 1:size(clin_df)[1] ## make so that no "missing" in TOP_ID or TOM_ID
    if ismissing(clin_df.TOP_ID[i])
        clin_df.TOP_ID[i] = "NA_"*string(i)
    end

    if ismissing(clin_df.TOM_ID[i])
        clin_df.TOM_ID[i] = "NA_"*string(i)
    end
end

### load and join metabolomics
linker_fh = project_path*mesa_data_sub_path*mesa_linker_fn
linker_df = CSV.read(linker_fh,DataFrame)

# HILIC LOAD AND TRANSPOSE
met_hil_fh = project_path*mesa_data_sub_path*mesa_hil_fn
met_hil_df = CSV.read(met_hil_fh,DataFrame,header = 5, skipto =6)
hil_names = names(met_hil_df)
met_hil_df.name = met_hil_df.Metabolite
hil_missing_name = ismissing.(met_hil_df.name)
met_hil_df.name[hil_missing_name] = met_hil_df.Compound[hil_missing_name]
met_hil_df_l = stack(met_hil_df, hil_names)
met_hil_df_t = unstack(met_hil_df_l, :variable, :name, :value)
met_hil_df_t = met_hil_df_t[7:end,:]
rename!(met_hil_df_t, 1 => :TOM_ID)
met_hil_df_t2 = passmissing(convert).(Float64,met_hil_df_t[:,2:end])
met_hil_df_t2.TOM_ID = met_hil_df_t.TOM_ID
met_hil_df_t = met_hil_df_t2

# AMIDE
met_ami_fh = project_path*mesa_data_sub_path*mesa_ami_fn
met_ami_df = CSV.read(met_ami_fh,DataFrame,footerskip = 2) ## ignoring two lines of notes at end of file

# JOIN
met_join_1 = outerjoin(met_hil_df_t, met_ami_df, on = :TOM_ID)
met_join = outerjoin(linker_df,met_join_1, on = :TOM_ID, makeunique = true)
filter!(row -> !ismissing(row.Exam), met_join)

### load proteomics
prot_fh = project_path*mesa_data_sub_path*mesa_prot_fn
prot_df = CSV.read(prot_fh,DataFrame,header = 8, skipto = 9)

clin_prot_join = outerjoin(clin_df,prot_df,on = :TOP_ID, makeunique = true)

mesa_join = outerjoin(clin_prot_join,met_join, on = :TOM_ID, makeunique = true)


### create arrays of variable names
met_ami_names = names(met_ami_df)[6:end]
met_hil_names = names(met_hil_df_t)[1:end-1]
met_names = vcat(met_ami_names,met_hil_names)

prot_names = names(prot_df)[34:end]

### Rename age and sex variables (sex.y and age.y are identical to the .x variables)
rename!(mesa_join, Symbol("gender1") => "sex")
rename!(mesa_join, Symbol("age1c") => "age")

### Rename BMI and eGFR variables
rename!(mesa_join, Symbol("bmi1c") => "bmi")
rename!(mesa_join, Symbol("cepgfr1t") => "egfr")

## Filter for Exam 1
filter!(row -> row.Exam == 1, mesa_join)

# remove names with significant missing data (>50%)
remove_hil_names = Set()
l = size(mesa_join)[1]
for met in met_hil_names
    if sum(ismissing.(mesa_join[:,Symbol(met)]))/l > 0.75
        push!(remove_hil_names,met)
    end
end
met_hil_names = collect(setdiff(Set(met_hil_names),remove_hil_names))

## Save and Reload DF for proper typing
CSV.write("$(project_path)$(mesa_data_sub_path)/temp_mesa_join.csv",mesa_join)
mesa_join = CSV.read("$(project_path)$(mesa_data_sub_path)/temp_mesa_join.csv",DataFrame)



## Perform Adjustments
mesa_join[:,:batch] .= 1
adj_vars!(Symbol.(prot_names),mesa_join,0,1,:batch,:age)
adj_vars!(prot_names,mesa_join,0,1,999,:sex)

adj_vars!(Symbol.(met_ami_names),mesa_join,0,1,999,:age)
adj_vars!(met_ami_names,mesa_join,0,1,999,:sex)

adj_vars!(Symbol.(met_hil_names),mesa_join,0,1,:batch,:age)
adj_vars!(met_hil_names,mesa_join,0,1,999,:sex)

CSV.write("$(project_path)$(mesa_data_sub_path)mesa_join_age_sex_norm.csv",mesa_join)

adj_vars!(prot_names,mesa_join,0,1,999,:bmi)
adj_vars!(met_ami_names,mesa_join,0,1,999,:bmi)
adj_vars!(met_hil_names,mesa_join,0,1,999,:bmi)

CSV.write("$(project_path)$(mesa_data_sub_path)mesa_join_age_sex_bmi_norm.csv",mesa_join)


adj_vars!(prot_names,mesa_join,0,1,999,:egfr)
adj_vars!(met_ami_names,mesa_join,0,1,999,:egfr)
adj_vars!(met_hil_names,mesa_join,0,1,999,:egfr)

CSV.write("$(project_path)$(mesa_data_sub_path)mesa_join_age_sex_bmi_egfr_norm.csv",mesa_join)


### write names files
prot_names_file = "$(project_path)$(mesa_data_sub_path)mesa_prot_names.csv"
f = open(prot_names_file, "w")

for prot in prot_names
	println(f, "\""*prot*"\"")
end
close(f)

met_names_file = "$(project_path)$(mesa_data_sub_path)mesa_met_names.csv"
f = open(met_names_file, "w")

for met in met_names
	println(f, "\""*met*"\"")
end
close(f)

ami_names_file = "$(project_path)$(mesa_data_sub_path)mesa_ami_names.csv"
f = open(ami_names_file, "w")

for ami in met_ami_names
	println(f, "\""*ami*"\"")
end
close(f)

hil_names_file = "$(project_path)$(mesa_data_sub_path)mesa_hil_names.csv"
f = open(hil_names_file, "w")

for hil in met_hil_names
	println(f, "\""*hil*"\"")
end
close(f)