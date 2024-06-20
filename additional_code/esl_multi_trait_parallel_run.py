import argparse
import os
import subprocess
from pathlib import Path
import shutil

def run_esl_analysis(input_dir, input_alignments_path, esl_directory, species_groups_file, job_name, rep_dir_name, param_dir_name, outgroup_species=None):
    user = os.environ.get('USER')
    
    species_combination_id = Path(species_groups_file).stem

    unique_tmp_dir = f"/local_scratch/tmp/{user}/{job_name}/{species_combination_id}/{rep_dir_name}/{param_dir_name}"
    esl_inputs_outputs_dir = f"{unique_tmp_dir}/esl_inputs_outputs"
    canceled_alignments_dir = f"{unique_tmp_dir}/canceled_alignments"

    os.makedirs(esl_inputs_outputs_dir, exist_ok=True)
    os.makedirs(canceled_alignments_dir, exist_ok=True)

    shutil.copytree(input_alignments_path, f"{unique_tmp_dir}/input_alignments", dirs_exist_ok=True)

    for simulated_file in Path(input_dir).glob("*_NT_simulated.fas"):
        base_name = simulated_file.stem.replace("_NT_simulated", "_AA")
        original_file_path = f"{unique_tmp_dir}/input_alignments/{base_name}.fas"
        try:
            os.remove(original_file_path)
        except FileNotFoundError:
            pass
        shutil.copy(simulated_file, f"{unique_tmp_dir}/input_alignments/")

    deletion_canceler_cmd = [
        "python", f"{esl_directory}/deletion_canceler.py",
        "--alignments_dir", f"{unique_tmp_dir}/input_alignments",
        "--canceled_alignments_dir", canceled_alignments_dir,
        "--species_groups_file", species_groups_file,
        "--cancel_only_partner",
        "--min_pairs", "3"
    ]
    
    if outgroup_species:
        deletion_canceler_cmd.extend(["--outgroup_species", outgroup_species])
    
    subprocess.run(deletion_canceler_cmd, check=True)

    subprocess.run([
        "python", f"{esl_directory}/esl_multimatrix.py",
        "--esl_inputs_outputs_dir", esl_inputs_outputs_dir,
        "--species_groups_file", species_groups_file,
        "--group_penalty_type", "median",
        "--use_logspace",
        "--output_file_base_name", f"{job_name}_{species_combination_id}_{rep_dir_name}_{param_dir_name}",
        "--no_pred_output",
        "--output_dir", f"{unique_tmp_dir}",
        "--use_existing_alignments",
        "--canceled_alignments_dir", f"{canceled_alignments_dir}",
        "--num_log_points", "4",
        "--delete_preprocess"
    ], check=True)

    output_csv = f"{unique_tmp_dir}/{job_name}_{species_combination_id}_{rep_dir_name}_{param_dir_name}_gene_ranks.csv"
    shutil.copy(output_csv, f"{Path(input_dir).parent}/{job_name}_{species_combination_id}_{rep_dir_name}_{param_dir_name}_gene_ranks.csv")

    shutil.rmtree(unique_tmp_dir)

def main():
    parser = argparse.ArgumentParser(description="Run ESL analysis for a single parameter directory of simulated alignments.")
    parser.add_argument("--input_dir", required=True, help="Directory containing the unique simulated alignments for the run.")
    parser.add_argument("--input_alignments_path", required=True, help="Path to the directory with input alignments before gap cancellation.")
    parser.add_argument("--esl_directory", required=True, help="Directory where the ESL scripts are located.")
    parser.add_argument("--species_groups_file", required=True, help="File specifying the species groups, used to derive unique identifier for temporary directories.")
    parser.add_argument("--job_name", required=True, help="Job name for the HPC queue.")
    parser.add_argument("--rep_dir_name", required=True, help="Name of the replicate directory.")
    parser.add_argument("--param_dir_name", required=True, help="Name of the parameter directory.")
    parser.add_argument("--outgroup_species", help="Optional outgroup species to be used in the analysis.")
    
    args = parser.parse_args()

    run_esl_analysis(args.input_dir, args.input_alignments_path, args.esl_directory, args.species_groups_file, args.job_name, args.rep_dir_name, args.param_dir_name, args.outgroup_species)

if __name__ == '__main__':
    main()

