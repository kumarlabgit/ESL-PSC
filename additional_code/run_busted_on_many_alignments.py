import os
import subprocess
import argparse
import json
from Bio import AlignIO

def alignment_length_check(alignment_path):
    """
    Checks if the length of sequences in an alignment exceeds a threshold.
    """
    alignment = AlignIO.read(alignment_path, "fasta")
    if len(alignment[0].seq) > 6000:
        return False
    return True

def generate_commands(alignment_dir, list_path, tree_path, commands_file_path):
    """
    Generate a command for each alignment file and write to a commands file.
    """
    script_path = "/home/user/work/scripts/busted_parallel_run.py"
    commands = []
    for filename in os.listdir(alignment_dir):
        if filename.endswith('.fas'):
            alignment_path = os.path.join(alignment_dir, filename)
            if alignment_length_check(alignment_path):
                command = f"timeout 5h python {script_path} {alignment_path} {tree_path} {list_path}"
                commands.append(command)
    
    with open(commands_file_path, 'w') as file:
        file.write("\n".join(commands))

def create_and_submit_job_script(commands_file_path, job_name, email, nodes, walltime):
    job_script_content = f"""#!/bin/sh
#PBS -N {job_name}
#PBS -l walltime={walltime}
#PBS -l nodes={nodes}:ppn=28
#PBS -q normal
#PBS -m bae
#PBS -M {email}

module load python/3.8.5
source /home/user/env_name/bin/activate

echo "Starting analysis with job name {job_name}."

timeout 13h torque-launch {commands_file_path}

echo "BUSTED commands finished. Starting to aggregate data and cleanup."

results_file="/home/user/work/busted_results_{job_name}.csv"
echo "alignment_name,p-value" > "$results_file"

for h in $(uniq $PBS_NODEFILE); do
    echo "Processing node: $h"
    ssh $h "cat /local_scratch/tmp/$USER/*_pvalue.txt" >> "$results_file"
    ssh $h "rm -rf /local_scratch/tmp/$USER" 
done

echo "Results aggregation completed. Results written to $results_file"
"""

    script_filename = f"{job_name}_job_script.sh"
    with open(script_filename, 'w') as file:
        file.write(job_script_content)
    
    subprocess.run(["qsub", script_filename])

def main(alignment_dir, list_path, tree_path, commands_file_name, job_name, email, nodes, walltime):
    commands_file_path = os.path.join(os.getcwd(), commands_file_name)
    generate_commands(alignment_dir, list_path, tree_path, commands_file_path)
    create_and_submit_job_script(commands_file_path, job_name, email, nodes, walltime)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and submit a job for running BUSTED analyses in parallel.")
    parser.add_argument("alignment_dir", help="Directory containing alignment files.")
    parser.add_argument("list_path", help="Path to the list of species.")
    parser.add_argument("tree_path", help="Path to the tree file.")
    parser.add_argument("--job_name", default="BUSTED_Analysis", help="Name of the job.")
    parser.add_argument("--email", required=True, help="Email for job notifications.")
    parser.add_argument("--nodes", default=1, help="Number of nodes to use.")
    parser.add_argument("--walltime", default="24:00:00", help="Walltime for the job.")
    
    args = parser.parse_args()
    main(args.alignment_dir, args.list_path, args.tree_path, args.job_name + "_commands.txt", args.job_name, args.email, args.nodes, args.walltime)

