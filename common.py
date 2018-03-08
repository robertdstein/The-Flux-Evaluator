from settings import root_path, source_path

root_path = root_path
source_path = source_path

tfe_path = root_path + "The-Flux-Evaluator__Data/"

input_path = tfe_path + "Input/"
storage_path = tfe_path + "Storage/"
output_path = tfe_path + "Output/"
log_path = tfe_path + "Logs/"

results_path = storage_path + "results"
merge_path = storage_path + "merged/results"
cat_path = input_path + "Catalogues/"
plot_path = output_path + "plots/"

coenders_7year_sens_path = source_path + "source/coenders_7year_sens.npy"

pickle_results_dir = output_path + "pickle_results/"

single_TDE_source_dir = cat_path + "Individual_TDEs/"
tde_pickle = pickle_results_dir + "Individual_TDEs.pkl"

tfe_sens_path = output_path + "PS_Sensitivities/"
tfe_k_sens_path = output_path + "PS_k_Sensitivities/"
tfe_seasons = ["IC40", "IC59", "IC79", "IC86_1", "IC86_234"]

new_data_dir = input_path + "PS_Data_Subsamples/"
e_decades_dir = new_data_dir + "Energy_Decades/"
