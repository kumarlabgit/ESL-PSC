import argparse
import os
from pathlib import Path

def generate_ccs_commands(top_level_dir, species_name, matrix_file):
    top_level_path = Path(top_level_dir)
    command_list = []

    for rep_dir in top_level_path.iterdir():
        if rep_dir.is_dir():
            for param_dir in rep_dir.iterdir():
                if param_dir.is_dir():
                    alignments_dir = str(param_dir)
                    command = f"python ccs_method.py {alignments_dir} {species_name} {matrix_file}"
                    command_list.append(command)
    return command_list

def submit_job_script(job_name, email, commands_file):
    job_script_content = f"""#!/bin/sh
#PBS -N {job_name}
#PBS -l walltime=1:00:00
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -m bae
#PBS -M {email}

cd $PBS_O_WORKDIR

module load python/3.8.5
source ~/env_name/bin/activate

torque-launch {commands_file}
"""

    script_filename = f"{job_name}_job_script.sh"
    with open(script_filename, "w") as file:
        file.write(job_script_content)
    os.system(f"qsub {script_filename}")

def main(args):
    commands = generate_ccs_commands(args.top_level_dir, args.species_name, args.matrix_file)
    
    commands_file = os.path.join(args.top_level_dir, 'ccs_commands.txt')
    with open(commands_file, 'w') as f:
        f.write("\n".join(commands))
    
    submit_job_script(args.job_name, args.email, commands_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate and submit jobs for CCS method.")
    parser.add_argument("top_level_dir", help="Top level directory containing the output from simulations.")
    parser.add_argument("species_name", help="Name of the outgroup species for CCS analysis.")
    parser.add_argument("matrix_file", help="Path to the matrix file (list of species in paired order).")
    parser.add_argument("--job_name", default="ccs_analysis", help="Job name for the HPC queue.")
    parser.add_argument("--email", required=True, help="Email for job notifications.")
    args = parser.parse_args()
    main(args)

