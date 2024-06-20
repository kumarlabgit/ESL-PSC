import argparse
import os
import subprocess
from pathlib import Path

def calculate_resources(num_alignments, time_per_alignment, cores_per_node=28, max_nodes=12):
    total_time_minutes = num_alignments * time_per_alignment / 60
    max_walltime_hours = 47 + 59 / 60  

    total_cores_needed = num_alignments  
    nodes_needed = min(-(-total_cores_needed // cores_per_node), max_nodes)

    total_batches = -(-total_cores_needed // (nodes_needed * cores_per_node))
    total_time_minutes = total_batches * time_per_alignment / 60

    if total_time_minutes / 60 > max_walltime_hours:
        estimated_walltime = f"{int(total_time_minutes // 60):02d}:{int(total_time_minutes % 60):02d}:00"
        print(f"Warning: Calculated wall time {estimated_walltime} exceeds the maximum limit of 47:59:00.")
        raise ValueError("Calculated wall time exceeds the maximum limit of 47:59:00.")

    walltime = f"{int(total_time_minutes // 60):02d}:{int(total_time_minutes % 60):02d}:00"
    return nodes_needed, walltime


def generate_commands(sim_output_dir, trees_dir, foreground_file):
    sim_output_path = Path(sim_output_dir)
    tree_path = Path(trees_dir)
    command_list = []

    for alignment_file in sim_output_path.rglob('*_NT_simulated.fas'):
        gene_name = alignment_file.stem.replace('_NT_simulated', '')
        tree_file = tree_path / f"{gene_name}_AA_pruned_tree.nwk"

        param_combo_dir = alignment_file.parent.name  
        rep_dir = alignment_file.parents[1].name  
        unique_dir_name = f"{rep_dir}_{param_combo_dir}_{gene_name}"  
        temp_dir = f"$TMPDIR/{unique_dir_name}"
        temp_alignment_file = f"{temp_dir}/{alignment_file.name}"
        temp_output_file = f"{temp_dir}/[required_output_file]"  
        final_output_dir = alignment_file.parent  

        command = (
            f"mkdir -p {temp_dir} && "  
            f"cp '{alignment_file}' {temp_dir} && "  
            f"cd {temp_dir} && "  
            f"csubst analyze --alignment_file {temp_alignment_file} "
            f"--rooted_tree_file '{tree_file}' "
            f"--foreground '{foreground_file}' "
            f"--exhaustive_until 1 --max_arity 10 --iqtree_exe iqtree2 && "
            f"mv {temp_dir}/csubst_cb_2.tsv '{final_output_dir}/{gene_name}_csubst_cb_2.tsv'"  
        )

        command_list.append(command)

    return command_list


def submit_job_script(job_name, email, walltime, nodes_needed, commands_file):
    job_script_content = f"""#!/bin/sh
#PBS -l walltime={walltime}
#PBS -N {job_name}
#PBS -q normal
#PBS -l nodes={nodes_needed}:ppn=28
#PBS -m bae
#PBS -M {email}
#PBS
cd $PBS_O_WORKDIR

module load python/3.8.5
source ~/env_name/bin/activate

torque-launch {commands_file}
"""

    script_filename = f"{job_name}_job_script.sh"
    with open(script_filename, "w") as file:
        file.write(job_script_content)

    subprocess.run(["qsub", script_filename])

def main(args):
    commands = generate_commands(args.sim_output_dir, args.trees_dir, args.foreground_file)
    
    commands_file = os.path.join(args.sim_output_dir, 'csubst_analyze_commands.txt')
    with open(commands_file, 'w') as f:
        f.write("\n".join(commands))
    
    nodes_needed, walltime = calculate_resources(len(commands), args.time_per_alignment)
    
    submit_job_script(args.job_name, args.email, walltime, nodes_needed, commands_file)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate and submit jobs for csubst analyze.")
    parser.add_argument("sim_output_dir", help="Directory containing the output from simulations.")
    parser.add_argument("trees_dir", help="Directory containing the tree files.")
    parser.add_argument("foreground_file", help="File specifying the foreground species.")
    parser.add_argument("--job_name", default="csubst_analyze", help="Job name for the HPC queue.")
    parser.add_argument("--email", required=True, help="Email for job notifications.")
    parser.add_argument("--time_per_alignment", type=int, default=60, help="Estimated time per alignment in seconds.")
    args = parser.parse_args()
    main(args)

