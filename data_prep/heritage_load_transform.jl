using DataFrames,CSV,StatsBase,RCall
try
    include("../analysis/omics_tools.jl")
catch
    include("./analysis/omics_tools.jl")
end


### Load data paths from data_prep.conf config file
conf = "./data_prep/data_prep.conf"
project_path = read_conf(conf,"project_path")
her_data_sub_path = read_conf(conf,"her_data_sub_path")
her_clin_fn = read_conf(conf,"her_clin_fn")
her_hil_fn = read_conf(conf,"her_hil_fn")
her_ami_ca_fn = read_conf(conf,"her_ami_ca_fn")
her_ami_aa_fn = read_conf(conf,"her_ami_aa_fn")
her_met_runlist_fn = read_conf(conf,"her_met_runlist_fn")
her_prot_fn = read_conf(conf,"her_prot_fn")

### join all heritage files into single table on PID
# load clinical
clin_fh = project_path*her_data_sub_path*her_clin_fn
clin_df = CSV.read(clin_fh,DataFrame)

### load met_runlist
her_met_runlist_fh = project_path*her_data_sub_path*her_met_runlist_fn
met_runlist = CSV.read(her_met_runlist_fh,DataFrame)
met_runlist_dict = Dict()
for i in 1:size(met_runlist)[1]
    run_number = string(met_runlist[i,Symbol("Run #")])
    pt_id = met_runlist[i,Symbol("pt ID")]
    timepoint = rstrip(met_runlist[i,Symbol("timepoint ")])
    met_runlist_dict[run_number] = [pt_id,timepoint]
end
rename!(met_runlist, 2 => :PID)
hil_batch_1_set = Set(string.(met_runlist.PID))
met_runlist.run = string.(met_runlist[:,Symbol("Run #")])

### load metabolites
# hilic load and transpose
met_hil_fh = project_path*her_data_sub_path*her_hil_fn
met_hil_df = CSV.read(met_hil_fh,DataFrame, header = 9, skipto = 10)
hil_names = names(met_hil_df)
met_hil_df.name = met_hil_df.Metabolite
hil_missing_name = ismissing.(met_hil_df.name)
met_hil_df.name[hil_missing_name] = met_hil_df.CompoundID_Heritage_1[hil_missing_name]
met_hil_df_l = stack(met_hil_df, hil_names)
met_hil_df_t = unstack(met_hil_df_l, :variable, :name, :value, allowduplicates = true)
met_hil_df_t = met_hil_df_t[10:end,:]
# replace all run list numbers with pt_id_timepoint
rename!(met_hil_df_t, 1 => :sample)
filter!(row -> !in(split(row.sample,"_")[1], hil_batch_1_set), met_hil_df_t)
met_hil_df_t.sample = rstrip.(met_hil_df_t.sample)
met_hil_df_t[!,:pt_id] .= ""
met_hil_df_t[!,:timepoint] .= ""
for i in 1:size(met_hil_df_t)[1]
    t_sample = met_hil_df_t.sample[i]
    if haskey(met_runlist_dict,t_sample)
        met_hil_df_t.sample[i] = join(met_runlist_dict[t_sample],"_")
    end
    sample_arr = split(met_hil_df_t.sample[i],"_")
    if length(sample_arr) == 2
        met_hil_df_t.pt_id[i] = sample_arr[1]
        met_hil_df_t.timepoint[i] = sample_arr[2]
    end
end

### need to remove second hilic datapoint for duplicate pateints run again in batch 2

# create only pt id column, and only timepoint column

# amide
met_ami_ca_fh = project_path*her_data_sub_path*her_ami_ca_fn
met_ami_ca_df = CSV.read(met_ami_ca_fh,DataFrame)
met_ami_ca_df[!,:batch] .= 1
for name in names(met_ami_ca_df)
    new_name = rstrip(replace(name," Results" => ""))
    rename!(met_ami_ca_df,Symbol(name) => Symbol(new_name))
end

met_ami_aa_fh = project_path*her_data_sub_path*her_ami_aa_fn
met_ami_aa_df = CSV.read(met_ami_aa_fh,DataFrame)
met_ami_aa_df[!,:batch] .= 2

#ca_names = Set(names(met_ami_ca_df))
#aa_names = Set(names(met_ami_aa_df))
rename!(met_ami_ca_df,Symbol("Succinic acid") => Symbol("Succinic acid/Methylmalonic acid"))
rename!(met_ami_ca_df,Symbol("pt ID") => Symbol("pt_id"))
rename!(met_ami_aa_df,Symbol("PID") => Symbol("pt_id"))
rename!(met_ami_aa_df,Symbol("Timept") => Symbol("timepoint"))

met_ami_ca_df[!,:sample] = rstrip.(string.(met_ami_ca_df.pt_id) .* "_" .* met_ami_ca_df.timepoint)

met_ami_aa_df[!,:sample] = string.(met_ami_aa_df.pt_id) .* "_" .* met_ami_aa_df.timepoint
met_ami_aa_df.sample[ismissing.(met_ami_aa_df.sample)] = met_ami_aa_df[ismissing.(met_ami_aa_df.pt_id),Symbol("Run #")]
met_ami_aa_df.sample = rstrip.(met_ami_aa_df.sample)

# remove ami bridging samples
ca_set = Set(met_ami_ca_df.pt_id)
filter!(row -> !in(row.pt_id,ca_set),met_ami_aa_df)

ami_comb = vcat(met_ami_ca_df,met_ami_aa_df,cols=:union)
filter(row -> endswith(row.sample,"Pre"),ami_comb) # 674
filter(row -> endswith(row.sample,"Pre"),met_hil_df_t) # 673

met_join = outerjoin(met_hil_df_t, ami_comb, on = :sample, makeunique = true)


#run_met_join = outerjoin(met_runlist,met_join,on = :run, makeunique = true)
filter!(row -> !ismissing(row.pt_id), met_join) # remove control samples without PIDs
filter!(row -> row.timepoint == "Pre", met_join) # keep only "Pre" samples ## 673

clin_df.pt_id = string.(clin_df.PID)
clin_met_join = rightjoin(clin_df,met_join,on = :pt_id)### APPEARS LIKE THERE ARE 15 PIDs WITHOUT CLINICAL INFORMATION

clin_met_join = innerjoin(clin_df,met_join,on = :pt_id)### do INNER join to remove patients without clinical information, 658

clin_pt_id_set = Set(string.(clin_df.PID))
met_pt_id_set = Set(met_join.pt_id)

# load proteins
her_prot_fh = project_path*her_data_sub_path*her_prot_fn
prot_df = CSV.read(her_prot_fh,DataFrame)

delete!(prot_df,765)
prot_df.pt_id = string.(prot_df.PID)
her_join = leftjoin(clin_met_join,prot_df,on = :pt_id, makeunique = true) # 658

## pt_id sets

prot_pt_id_set = Set(string.(collect(prot_df.pt_id)))

intersect(clin_pt_id_set,prot_pt_id_set)

### create arrays of variable names
met_ami_names = names(ami_comb)[4:end]
met_ami_names = setdiff(met_ami_names,Set(["sample","Samples","location","batch","Box"]))

met_hil_names = names(met_hil_df_t)[2:end-2]
met_names = vcat(met_ami_names,met_hil_names)
# remove names with significant missing data (>50%)
remove_hil_names = Set()
l = size(met_hil_df_t)[1]
for met in met_hil_names
    if sum(ismissing.(met_hil_df_t[:,Symbol(met)]))/l > 0.50
        push!(remove_hil_names,met)
    end
end
met_hil_names = collect(setdiff(Set(met_hil_names),remove_hil_names))

prot_names = names(prot_df)[3:end-1]
prot_names = prot_names[occursin.(r"^B_",prot_names)] # only include baseline


### Rename age and sex variables (sex.y and age.y are identical to the .x variables)
rename!(her_join, Symbol("SEX") => "sex")
rename!(her_join, Symbol("AGE") => "age")

### Rename BMI and eGFR variables
rename!(her_join, Symbol("B_BMI") => "bmi")

### Determine patient rows with complete protein data, modify jhs_df to only include those rows --->> change to check for complete cases for each pairwise comparison, not feasible to have complete cases across the board when using 5000 plex platform

req_vars = vcat([:age,:sex,:bmi],Symbol.(prot_names))


#her_join[:,:batch] .= 1

### manual data file edits
# CA - remove single negative value for N-Acetyl-Glutamic-acid, 2 negative values for UDP-glucose / -galactose
# AA - remove single negative value for oxalic acid

adj_vars!(Symbol.(prot_names),her_join,0,1,:batch,:age) ## assume proteins were done in same ca/aa batches as the amide metabolites
adj_vars!(prot_names,her_join,0,1,999,:sex)

#her_join[!,Symbol("N-Acetyl-L-Glutamic acid")]
#her_join[!,Symbol("N-Acetyl-L-Aspartic acid")]

adj_vars!(Symbol.(met_ami_names),her_join,0,1,:batch,:age)
adj_vars!(met_ami_names,her_join,0,1,999,:sex)

her_join[!,met_hil_names] = passmissing(convert).(Float64,her_join[:,met_hil_names])
#CSV.write("$(project_path)$(her_data_sub_path)temp_her_join.csv",her_join)
#her_join = CSV.read("$(project_path)$(her_data_sub_path)temp_her_join.csv",DataFrame)    

adj_vars!(Symbol.(met_hil_names),her_join,0,1,:batch,:age)
adj_vars!(met_hil_names,her_join,0,1,999,:sex)

CSV.write("$(project_path)$(her_data_sub_path)her_join_age_sex_norm.csv",her_join)

adj_vars!(prot_names,her_join,0,1,999,:bmi)
adj_vars!(met_ami_names,her_join,0,1,999,:bmi)
adj_vars!(met_hil_names,her_join,0,1,999,:bmi)

CSV.write("$(project_path)$(her_data_sub_path)her_join_age_sex_bmi_norm.csv",her_join)


### write names files
prot_names_file = "$(project_path)$(her_data_sub_path)her_prot_names.csv"
f = open(prot_names_file, "w")

for prot in prot_names
	println(f, "\""*prot*"\"")
end
close(f)

met_names_file = "$(project_path)$(her_data_sub_path)her_met_names.csv"
f = open(met_names_file, "w")

for met in met_names
	println(f, "\""*met*"\"")
end
close(f)

ami_names_file = "$(project_path)$(her_data_sub_path)her_ami_names.csv"
f = open(ami_names_file, "w")

for ami in met_ami_names
	println(f, "\""*ami*"\"")
end
close(f)

hil_names_file = "$(project_path)$(her_data_sub_path)her_hil_names.csv"
f = open(hil_names_file, "w")

for hil in met_hil_names
	println(f, "\""*hil*"\"")
end
close(f)