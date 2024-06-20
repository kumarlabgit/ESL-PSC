import argparse
import os
import subprocess
from pathlib import Path
import shutil

def run_esl_analysis(input_dir, background_alignments_path, esl_directory, species_groups_file, job_name, rep_dir_name, param_dir_name):
    user = os.environ.get('USER')
    
    unique_tmp_dir = f"/local_scratch/tmp/{user}/{job_name}/{rep_dir_name}/{param_dir_name}"
    simulated_alignments_dir = f"{unique_tmp_dir}/simulated_alignments"
    esl_inputs_outputs_dir = f"{unique_tmp_dir}/esl_inputs_outputs"
    canceled_simulated_dir = f"{unique_tmp_dir}/gap-canceled_simulated"

    os.makedirs(esl_inputs_outputs_dir, exist_ok=True)
    os.makedirs(canceled_simulated_dir, exist_ok=True)

    shutil.copytree(background_alignments_path, f"{unique_tmp_dir}/background_alignments", dirs_exist_ok=True)
    shutil.copytree(input_dir, simulated_alignments_dir, dirs_exist_ok=True)

    subprocess.run([
        "python", f"{esl_directory}/deletion_canceler.py",
        "--alignments_dir", simulated_alignments_dir,
        "--canceled_alignments_dir", canceled_simulated_dir,
        "--species_groups_file", species_groups_file
    ], check=True)

    for i in range(16):
        combo_dir = f"combo_{i}-alignments"
        for file_path in Path(canceled_simulated_dir).glob(f"{combo_dir}/*_NT_simulated.fas"):
            base_name = file_path.stem.replace("_NT_simulated", "_AA")
            original_file_path = f"{unique_tmp_dir}/background_alignments/{combo_dir}/{base_name}.fas"
            
            try:
                os.remove(original_file_path)
            except FileNotFoundError:
                pass  
            
            shutil.copy(file_path, f"{unique_tmp_dir}/background_alignments/{combo_dir}/")

    subprocess.run([
        "python", f"{esl_directory}/esl_multimatrix.py",
        "--esl_inputs_outputs_dir", esl_inputs_outputs_dir,
        "--species_groups_file", species_groups_file,
        "--group_penalty_type", "median",
        "--use_logspace",
        "--output_file_base_name", f"{job_name}_{rep_dir_name}_{param_dir_name}",
        "--no_pred_output",
        "--output_dir", f"{unique_tmp_dir}",
        "--num_log_points", "4",
        "--use_existing_alignments",
        "--canceled_alignments_dir", f"{unique_tmp_dir}/background_alignments",
        "--delete_preprocess"
    ], check=True)

    output_csv = f"{unique_tmp_dir}/{job_name}_{rep_dir_name}_{param_dir_name}_gene_ranks.csv"
    shutil.copy(output_csv, f"{Path(input_dir).parent}/{job_name}_{rep_dir_name}_{param_dir_name}_gene_ranks.csv")

    shutil.rmtree(unique_tmp_dir)

def main():
    parser = argparse.ArgumentParser(description="Run ESL analysis for a single parameter directory of simulated alignments.")
    parser.add_argument("--input_dir", required=True, help="Directory containing the unique simulated alignments for the run.")
    parser.add_argument("--background_alignments_path", required=True, help="Path to the directory with background alignments.")
    parser.add_argument("--esl_directory", required=True, help="Directory where the ESL scripts are located.")
    parser.add_argument("--species_groups_file", required=True, help="File specifying the species groups.")
    parser.add_argument("--job_name", required=True, help="Job name for the HPC queue.")
    parser.add_argument("--rep_dir_name", required=True, help="Name of the replicate directory.")
    parser.add_argument("--param_dir_name", required=True, help="Name of the parameter directory.")
    
    args = parser.parse_args()

    run_esl_analysis(args.input_dir, args.background_alignments_path, args.esl_directory, args.species_groups_file, args.job_name, args.rep_dir_name, args.param_dir_name)

if __name__ == '__main__':
    main()

