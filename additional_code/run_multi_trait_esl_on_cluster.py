import argparse
import subprocess
from pathlib import Path

def generate_esl_commands(species_groups_dir, input_dir, esl_directory, job_name, root_dir_name, outgroup_species=None):
    commands = []
    species_groups_files = {file.stem: file for file in Path(species_groups_dir).glob('*.txt')}

    for combo_dir in Path(input_dir).glob('*'):
        combo_name = combo_dir.name
        if combo_dir.is_dir() and combo_name in species_groups_files:
            species_groups_file = species_groups_files[combo_name]
            for rep_dir in combo_dir.glob('rep*'):
                if rep_dir.is_dir():
                    for param_dir in rep_dir.glob('*'):
                        if param_dir.is_dir():
                            command = (
                                f"python /home/user/work/scripts/esl_multi_trait_parallel_run.py "
                                f"--input_dir {param_dir} "
                                f"--input_alignments_path /local_scratch/tmp/$USER/{job_name}/{root_dir_name} "
                                f"--esl_directory {esl_directory} "
                                f"--species_groups_file {species_groups_file} "
                                f"--job_name {job_name} "
                                f"--rep_dir_name {rep_dir.name} "
                                f"--param_dir_name {param_dir.name}"
                            )
                            if outgroup_species:
                                command += f" --outgroup_species {outgroup_species}"
                            commands.append(command)
    
    commands_file_path = Path(input_dir).parent / f"{job_name}_esl_commands.txt"
    with open(commands_file_path, 'w') as f:
        for command in commands:
            f.write(command + "\n")
            
    return commands_file_path

def create_and_submit_job_script(commands_file_path, job_name, email, input_alignments_tar_gz, root_dir_name, nodes, walltime):
    job_script_content = f"""#!/bin/bash
#PBS -N {job_name}
#PBS -l walltime={walltime}
#PBS -l nodes={nodes}:ppn=16
#PBS -q medium
#PBS -m bae
#PBS -M {email}

module load python/3.8.5
source /home/user/env_name/bin/activate

for h in $(uniq < $PBS_NODEFILE); do
    ssh $h "mkdir -p /local_scratch/tmp/$USER/{job_name}"
    scp "{input_alignments_tar_gz}" $h:/local_scratch/tmp/$USER/{job_name}/
    ssh $h "tar -xzvf /local_scratch/tmp/$USER/{job_name}/$(basename "{input_alignments_tar_gz}") -C /local_scratch/tmp/$USER/{job_name}"
done

torque-launch {commands_file_path}

for h in $(uniq $PBS_NODEFILE); do
    ssh $h "rm -rf /local_scratch/tmp/$USER/{job_name}"
done
"""

    script_filename = f"{job_name}_job_script.sh"
    with open(script_filename, 'w') as file:
        file.write(job_script_content.format(
            job_name=job_name, 
            input_alignments_tar_gz=input_alignments_tar_gz, 
            commands_file_path=commands_file_path,
            walltime=walltime,
            nodes=nodes,
            email=email,
            root_dir_name=root_dir_name))

    subprocess.run(["qsub", script_filename])

def main():
    parser = argparse.ArgumentParser(description="Script to submit ESL analysis jobs on HPC for multiple species sets.")
    parser.add_argument("--species_groups_dir", required=True, help="Directory of species groups files corresponding to each species combination.")
    parser.add_argument("--input_dir", required=True, help="Top-level directory containing all the combination directories with simulated sets.")
    parser.add_argument("--input_alignments_tar_gz", required=True, help="Path to the .tar.gz file with the full input alignments.")
    parser.add_argument("--root_dir_name", required=True, help="Name of the root directory in the input alignments tarball.")
    parser.add_argument("--esl_directory", required=True, help="Directory where the ESL scripts are located.")
    parser.add_argument("--job_name", default="ESL_MultiTrait_Analysis", help="Job name for the HPC queue.")
    parser.add_argument("--email", required=True, help="Email for job notifications.")
    parser.add_argument("--nodes", default=8, type=int, help="Number of nodes to request.")
    parser.add_argument("--walltime", default="18:00:00", help="Wall time for the job (format: HH:MM:SS).")
    parser.add_argument("--outgroup_species", help="Optional outgroup species to be used in the analysis.")
    
    args = parser.parse_args()

    commands_file_path = generate_esl_commands(args.species_groups_dir, args.input_dir, args.esl_directory, args.job_name, args.root_dir_name, args.outgroup_species)
    create_and_submit_job_script(commands_file_path, args.job_name, args.email, args.input_alignments_tar_gz, args.root_dir_name, args.nodes, args.walltime)

if __name__ == '__main__':
    main()

