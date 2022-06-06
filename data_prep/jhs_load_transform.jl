using DataFrames, CSV, StatsBase, RCall
try
    include("../analysis/omics_tools.jl")
catch
    include("./analysis/omics_tools.jl")
end


### Load data paths from data_prep.conf config file
conf = "./data_prep/data_prep.conf"
project_path = read_conf(conf,"project_path")
jhs_data_sub_path = read_conf(conf,"jhs_data_sub_path")
jhs_data_fn = read_conf(conf,"jhs_data_fn")

jhs_path = project_path*jhs_data_sub_path*jhs_data_fn

### read in CSV file
jhs_df = CSV.read(jhs_path, DataFrame)

### remove variables with _impute and _binary name suffixes
for var in names(jhs_df)
    var_str = string(var)
    if occursin(r"_impute$|_binary$",var_str)
        select!(jhs_df,Not(var))
    end
end

## Remove proteins with only 399 values
select!(jhs_df,Not([:SL002077,:SL000311,:SL004857,:SL004791,:SL008309]))

### Determine protein column names ###
jhs_var_names = String.(names(jhs_df))
jhs_protein_names = jhs_var_names[(occursin.(r"^SL0",jhs_var_names))]
jhs_protein_sym = Symbol.(jhs_protein_names)

### Rename age and sex variables (sex.y and age.y are identical to the .x variables)
rename!(jhs_df, Symbol("sex.x") => "sex")
rename!(jhs_df, Symbol("age.x") => "age")

### Rename BMI and eGFR variables
rename!(jhs_df, Symbol("BMI.x") => "bmi")
rename!(jhs_df, Symbol("eGFR.x") => "egfr")

### Define Metabolite column names ###
jhs_met_names = jhs_var_names[2108:end]
# remove QI
jhs_met_names = jhs_met_names[.!occursin.(r"^QI",jhs_met_names)]
jhs_met_sym = Symbol.(jhs_met_names)

### Define all Variable names ###
prot_met_names = vcat(jhs_protein_names,jhs_met_names)
clin_vars = ["subjid","sex","age","bmi","egfr","prot_batch"]
all_names = vcat(clin_vars,prot_met_names)

### Select only all_names vars
jhs_df = jhs_df[:,Symbol.(all_names)]

### Filter rows without age
filter!(row -> .!ismissing.(row.age),jhs_df)


### Scale and Adjust Protein Variables for Age and Sex ###
#adj_vars!(jhs_protein_sym,jhs_df,0,1,:prot_batch,:age,:sex)

# Adjust prot and met for age/sex and log norm prot by batch
adj_vars!(jhs_protein_sym,jhs_df,0,1,:prot_batch,:age)
adj_vars!(jhs_protein_sym,jhs_df,0,1,999,:sex)
adj_vars!(jhs_met_sym,jhs_df,0,1,999,:age)
adj_vars!(jhs_met_sym,jhs_df,0,1,999,:sex)

CSV.write("$(project_path)$(jhs_data_sub_path)jhs_join_age_sex_norm.csv",jhs_df)

adj_vars!(jhs_protein_sym,jhs_df,0,1,999,:bmi)
adj_vars!(jhs_met_sym,jhs_df,0,1,999,:bmi)

CSV.write("$(project_path)$(jhs_data_sub_path)jhs_join_age_sex_bmi_norm.csv",jhs_df)

adj_vars!(jhs_protein_sym,jhs_df,0,1,999,:egfr)
adj_vars!(jhs_met_sym,jhs_df,0,1,999,:egfr)

CSV.write("$(project_path)$(jhs_data_sub_path)jhs_join_age_sex_bmi_egfr_norm.csv",jhs_df)

### write names files
prot_names_file = "$(project_path)$(jhs_data_sub_path)jhs_prot_names.csv"
f = open(prot_names_file, "w")

for prot in jhs_protein_names
	println(f, "\""*prot*"\"")
end
close(f)

met_names_file = "$(project_path)$(jhs_data_sub_path)jhs_met_names.csv"
f = open(met_names_file, "w")

for met in jhs_met_names
	println(f, "\""*met*"\"")
end
close(f)