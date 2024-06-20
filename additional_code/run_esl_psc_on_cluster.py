import argparse
import subprocess
from pathlib import Path

def generate_esl_commands(input_dir, esl_directory, species_groups_file, job_name, root_dir_name):
    commands = []
    background_alignments_path = f"/local_scratch/tmp/$USER/{job_name}/{root_dir_name}"
    
    for rep_dir in Path(input_dir).glob('rep*'):
        if rep_dir.is_dir():
            for param_dir in rep_dir.glob('*'):
                if param_dir.is_dir():
                    command = (
                        f"python /home/user/work/scripts/esl_parallel_run.py "
                        f"--input_dir {param_dir} "
                        f"--background_alignments_path {background_alignments_path} "
                        f"--esl_directory {esl_directory} "
                        f"--species_groups_file {species_groups_file} "
                        f"--job_name {job_name} "
                        f"--rep_dir_name {rep_dir.name} "
                        f"--param_dir_name {param_dir.name}"
                    )
                    commands.append(command)
    
    commands_file_path = Path(input_dir) / "esl_commands.txt"
    with open(commands_file_path, 'w') as f:
        for command in commands:
            f.write(command + "\n")
            
    return commands_file_path

def create_and_submit_job_script(commands_file_path, job_name, email, background_alignments_tar_gz, nodes, walltime):
    job_script_content = f"""#!/bin/bash
#PBS -N {job_name}
#PBS -l walltime={walltime}
#PBS -l nodes={nodes}:ppn=18
#PBS -q normal
#PBS -m bae
#PBS -M {email}

module load python/3.8.5
source ~/env_name/bin/activate

for h in $(uniq < $PBS_NODEFILE); do
    ssh $h "mkdir -p /local_scratch/tmp/$USER/{job_name}"
    scp "{background_alignments_tar_gz}" $h:/local_scratch/tmp/$USER/{job_name}/
    ssh $h "tar -xzvf /local_scratch/tmp/$USER/{job_name}/$(basename "${background_alignments_tar_gz}") -C /local_scratch/tmp/$USER/{job_name}"
done

torque-launch {commands_file_path}

for h in $(uniq $PBS_NODEFILE); do
    ssh $h "rm -rf /local_scratch/tmp/$USER/{job_name}"
done
"""

    script_filename = f"{job_name}.sh"
    with open(script_filename, 'w') as file:
        file.write(job_script_content.format(job_name=job_name, background_alignments_tar_gz=background_alignments_tar_gz, commands_file_path=commands_file_path))

    subprocess.run(["qsub", script_filename])

def main():
    parser = argparse.ArgumentParser(description="Script to submit ESL analysis jobs on HPC")
    parser.add_argument("--input_dir", required=True, help="Directory containing the unique simulated alignments for each run.")
    parser.add_argument("--background_alignments_tar_gz", required=True, help="Path to the .tar.gz file with background alignments.")
    parser.add_argument("--esl_directory", required=True, help="Directory where the ESL scripts are located.")
    parser.add_argument("--species_groups_file", required=True, help="File specifying the species groups.")
    parser.add_argument("--job_name", default="ESL_Analysis", help="Job name for the HPC queue.")
    parser.add_argument("--email", required=True, help="Email for job notifications.")
    parser.add_argument("--nodes", default=8, help="number of nodes to request")
    parser.add_argument("--walltime", default="18:00:00", help="Wall time for the job (format: HH:MM:SS). it takes around 1hr 45mins for 1 round of ESL-PSC on a medium node")    
    parser.add_argument("--root_dir_name", required=True, help="Name of the root directory in the background alignments tarball.")
    
    args = parser.parse_args()

    commands_file_path = generate_esl_commands(args.input_dir, args.esl_directory, args.species_groups_file, args.job_name, args.root_dir_name)
    create_and_submit_job_script(commands_file_path, args.job_name, args.email, args.background_alignments_tar_gz, args.nodes, args.walltime)

if __name__ == '__main__':
    main()
