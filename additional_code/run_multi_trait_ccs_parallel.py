import argparse
import os
from pathlib import Path

def generate_ccs_commands(species_groups_dir, top_level_dir, outgroup_species, ignore_missing_species=False, job_name="", background_alignments_root=""):
    command_list = []
    species_groups_files = {file.stem: file for file in Path(species_groups_dir).glob('*.txt')}
    background_output_dir = ""

    user = os.environ.get('USER')

    local_alignments_dir = f"/local_scratch/tmp/{user}/{job_name}/{background_alignments_root}" if background_alignments_root else None

    if background_alignments_root:
        background_output_dir = Path(top_level_dir).parent / f"{Path(species_groups_dir).name}_ccs_background"
        os.makedirs(background_output_dir, exist_ok=True)

    for combo_dir in Path(top_level_dir).glob('*'):
        combo_name = combo_dir.name
        if combo_dir.is_dir() and combo_name in species_groups_files:
            matrix_file = species_groups_files[combo_name]
            if background_alignments_root:
                output_path = background_output_dir / f"{combo_name}.csv"
                command = f"python ccs_method.py {local_alignments_dir} {' '.join(outgroup_species)} {matrix_file} --output_path {output_path}"
                if ignore_missing_species:
                    command += " --ignore_missing_species"
                command_list.append(command)
            for rep_dir in combo_dir.glob('rep*'):
                if rep_dir.is_dir():
                    for param_dir in rep_dir.glob('*'):
                        if param_dir.is_dir():
                            alignments_dir = str(param_dir)
                            command = f"python ccs_method.py {alignments_dir} {' '.join(outgroup_species)} {matrix_file}"
                            if ignore_missing_species:
                                command += " --ignore_missing_species"
                            command_list.append(command)
    return command_list, background_output_dir

def submit_job_script(job_name, email, commands_file, alignments_tar_gz_path="", background_output_dir=""):
    job_script_content = f"""#!/bin/sh
#PBS -N {job_name}
#PBS -l walltime=2:00:00
#PBS -q medium
#PBS -l nodes=1:ppn=16
#PBS -m bae
#PBS -M {email}
#PBS
cd $PBS_O_WORKDIR

module load python/3.8.5
source /home/user/env_name/bin/activate
"""

    if alignments_tar_gz_path:
        distribute_alignments = f"""
for h in $(uniq < $PBS_NODEFILE); do
    ssh $h "mkdir -p /local_scratch/tmp/$USER/{job_name}"
    scp "{alignments_tar_gz_path}" $h:/local_scratch/tmp/$USER/{job_name}/
    ssh $h "tar -xzvf /local_scratch/tmp/$USER/{job_name}/$(basename "{alignments_tar_gz_path}") -C /local_scratch/tmp/$USER/{job_name}"
done
"""
        job_script_content += distribute_alignments

    job_script_content += f"\ntorque-launch {commands_file}\n"

    if background_output_dir:
        job_script_content += f"\n

    script_filename = f"{job_name}_job_script.sh"
    with open(script_filename, "w") as file:
        file.write(job_script_content)
    os.system(f"qsub {script_filename}")

def main():
    parser = argparse.ArgumentParser(description="Generate and submit jobs for CCS method across multiple species sets.")
    parser.add_argument("species_groups_dir", help="Directory of species groups files for each species set.")
    parser.add_argument("top_level_dir", help="Top level directory containing the species set directories.")
    parser.add_argument("outgroup_species", nargs='+', help="List of the outgroup species IDs for CCS analysis.")
    parser.add_argument("--ignore_missing_species", action='store_true', help="Ignore missing species and gaps in the alignments.")
    parser.add_argument("--generate_background", help="Path to the tarball containing all the original alignments for background distribution.")
    parser.add_argument("--background_alignments_root", help="Name of the root directory within the background alignments tar.gz.")
    parser.add_argument("--job_name", default="ccs_multi_trait_analysis", help="Job name for the HPC queue.")
    parser.add_argument("--email", required=True, help="Email for job notifications.")
    args = parser.parse_args()

    commands, background_output_dir = generate_ccs_commands(args.species_groups_dir, args.top_level_dir, args.outgroup_species, args.ignore_missing_species, args.job_name, args.background_alignments_root)
    
    commands_file = os.path.join(args.top_level_dir, 'ccs_multi_trait_commands.txt')
    with open(commands_file, 'w') as f:
        f.write("\n".join(commands))
    
    submit_job_script(args.job_name, args.email, commands_file, args.generate_background if args.generate_background else "")

if __name__ == '__main__':
    main()



