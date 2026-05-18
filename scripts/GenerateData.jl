using PhysiCellModelManager
initializeModelManager() # this works if launching from the project directory, i.e. the directory containing the data and PhysiCell folders
# initializeModelManager("/Users/jeanette.johnson/PDACAntigenPresentationABMs/PhysiCell", "/Users/jeanette.johnson/PDACAntigenPresentationABMs/data") # use this if not calling this from the project directory

############ set up ############

config_folder = "0_template" # this folder is located at ./data/inputs/configs
custom_code_folder = "0_template" # this folder is located at ./data/inputs/custom_codes
rulesets_collection_folder = "0_template" # this folder is located at ./data/inputs/rulesets_collections; a rule has been added for the sake of the example
intracellular_folder = "" # optionally add this folder with intracellular.xml to ./data/inputs/ics/intracellulars

ic_cell_folder = "0_template" # this folder is located at ./data/inputs/ics/cells
ic_substrate_folder = "" # optionally add this folder with substrates.csv to ./data/inputs/ics/substrates
ic_ecm_folder = "" # optionally add this folder with ecms.csv to ./data/inputs/ics/ecms
ic_dc_folder = "" # optionally add this folder with dcs.csv to ./data/inputs/ics/dcs

# package them all together into a single object
inputs = InputFolders(config_folder, custom_code_folder;
                        rulesets_collection=rulesets_collection_folder,
                        intracellular=intracellular_folder,
                        ic_cell=ic_cell_folder,
                        ic_substrate=ic_substrate_folder,
                        ic_ecm=ic_ecm_folder,
                        ic_dc=ic_dc_folder)

############ make the simulations short ############

# We will set the default simulations to have a lower max time.
# This will serve as a reference for the following simulations.
xml_path = ["overall"; "max_time"]
value = 60.0
dv_max_time = DiscreteVariation(xml_path, value)
reference = createTrial(inputs, dv_max_time; n_replicates=0) # since we don't want to run this, set the n_replicates to 0

############ set up variables to control running simulations ############

# you can force the recompilation, but it is usually only necesary if you change core code
# if you change custom code, it is recommended you make a new custom codes folder in ./data/inputs/custom_codes...
# ...especially if the database already has simulations run with that custom code
force_recompile = false

# PhysiCellModelManager.jl records which simulations all use the same parameter vector...
# ...to reuse them (unless the user opts out)
use_previous = true # if true, will attempt to reuse simulations with the same parameters; otherwise run new simulations

# a monad refers to a single collection of identical simulations...
# except for randomness (could be do to the initial seed or stochasticity introduced by omp threading)
# n_replicates is the number of replicates to run for each parameter vector...
# ...PhysiCellModelManager.jl records which simulations all use the same parameter vector...
# ...and will attempt to reuse these (unless the user opts out)...
# ...so this parameter is the _min_ because there may already be many sims with the same parameters
n_replicates = 1

############ set up parameter variations ############

# assume you have the template project with "default" as a cell type...
# ...let's vary their cycle durations and apoptosis rates

# get the xml path to duration of phase 0 of the default cell type
# this is a list of strings in which each string is either...
# 	1) the name of a tag in the xml file OR
# 	2) the name of a tag along with the value of one attribute (name:attribute_name:attribute_value)
xml_path = PhysiCellModelManager.cyclePath("default", "phase_durations", "duration:index:0")
vals = [200.0, 300.0, 400.0] # choose 3 discrete values to vary the duration of phase 0
dv_phase_0_duration = DiscreteVariation(xml_path, vals)

# now do the same, but for the apoptosis rate
xml_path = PhysiCellModelManager.apoptosisPath("default", "death_rate")
vals = [4.31667e-05, 5.31667e-05, 6.31667e-05] # choose 3 discrete values to vary the apoptosis rate
dv_apoptosis_rate = DiscreteVariation(xml_path, vals)

# now combine them into a list:
discrete_variations = [dv_phase_0_duration, dv_apoptosis_rate]

############ run the sampling ############

# now create the sampling (varied parameter values) with these parameters
# we will give it a reference to the monad with the short max time
sampling = createTrial(reference, discrete_variations; n_replicates=n_replicates, use_previous=use_previous)

# at this point, we have only added the sampling to the database...
# ...along with the monads and simulations that make it up
# before running, we will set the number of parallel simulations to run.
# note: this will only be used when running locally, i.e., not on an HPC
# by default, PhysiCellModelManager.jl will run the simulations serially, i.e., 1 in "parallel".
# change this by calling:
setNumberOfParallelSims(4) # for example, to run 4 simulations in parallel

# you can change this default behavior on your machine by setting an environment variable...
# called PCMM_NUM_PARALLEL_SIMS
# this is read during `initializeModelManager`...
# meaning subsequent calls to `setNumberOfParallelSims` will overwrite the value
# A simple way to use this when running the script is to run in your shell:
# `PCMM_NUM_PARALLEL_SIMS=4 julia ./scripts/GenerateData.jl`

# now run the sampling
out = run(sampling; force_recompile=force_recompile)

# If you are running on an SLURM-based HPC, PhysiCellModelManager.jl will detect this and calls to `sbatch`...
# ...to parallelize the simulations, batching out each simulation to its own job.
