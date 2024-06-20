
import argparse
import subprocess
import os

def calculate_resources(num_replicates, num_align_per_rep, num_sites, num_scaling_factors, time_per_command):
    total_commands = num_replicates * num_align_per_rep * len(num_sites) * len(num_scaling_factors)
    cores_per_node = 28  
    max_nodes = 12  
    
    nodes_needed = min(-(-total_commands // cores_per_node), max_nodes)

    total_batches = -(-total_commands // (nodes_needed * cores_per_node))
    
    total_time_minutes = total_batches * time_per_command / 60

    max_walltime_hours = 47 + 59 / 60  

    if total_time_minutes / 60 > max_walltime_hours:
        raise ValueError("Calculated wall time exceeds the maximum limit of 47:59:00.")

    walltime = f"{int(total_time_minutes // 60):02d}:{int(total_time_minutes % 60):02d}:00"
    
    return nodes_needed, walltime


def submit_job_script(args):
    nodes, walltime = calculate_resources(args.num_replicates, args.num_align_per_rep, args.num_sites, args.scaling_factors, args.time_per_command)

    job_name = os.path.basename(args.output_dir.rstrip('/'))
    script_filename = f"{job_name}_job_script.sh"
    allowed_proteins_arg = f"--allowed_proteins_file {args.allowed_proteins_file}" if args.allowed_proteins_file else ""
    rep_lists_arg = f"--rep_lists {args.rep_lists}" if args.rep_lists else ""
    checkpoint_file = f"{job_name}_chkpt.txt"  

    job_script_content = f"""#!/bin/sh
#PBS -l walltime={walltime}
#PBS -N {job_name}
#PBS -A tutorial
#PBS -q normal
#PBS -l nodes={nodes}:ppn=28
#PBS -m bae
#PBS -M {args.email}
#PBS
cd $PBS_O_WORKDIR

module load python/3.8.5
source ~/env_name/bin/activate

python3 sim_preprocess_for_esl_benchmarks.py {args.alignments_dir} {','.join(map(str, args.num_sites))} {','.join(map(str, args.scaling_factors))} {args.trees_dir} {args.output_dir} {args.foreground_file} {args.convergent_aa} {args.num_align_per_rep} {args.num_replicates} {allowed_proteins_arg} {rep_lists_arg}

torque-launch -p {checkpoint_file} {os.path.join(args.output_dir, 'simulation_commands.txt')}
"""

    with open(script_filename, "w") as file:
        file.write(job_script_content)

    subprocess.run(["qsub", script_filename])

def parse_arguments():
    parser = argparse.ArgumentParser(description='HPC Job Submission Script')
    parser.add_argument('--alignments_dir', required=True, help='Path to alignments directory')
    parser.add_argument('--num_sites', type=lambda s: [int(item) for item in s.split(',')], required=True, help='List of numbers of sites (comma-separated)')
    parser.add_argument('--scaling_factors', type=lambda s: [float(item) for item in s.split(',')], required=True, help='List of scaling factors (comma-separated)')
    parser.add_argument('--trees_dir', required=True, help='Path to trees directory')
    parser.add_argument('--output_dir', required=True, help='Path to output directory')
    parser.add_argument('--foreground_file', required=True, help='Foreground file')
    parser.add_argument('--convergent_aa', required=True, help='Convergent amino acids')
    parser.add_argument('--num_align_per_rep', type=int, required=True, help='Number of alignments per replicate')
    parser.add_argument('--num_replicates', type=int, required=True, help='Number of replicates')
    parser.add_argument('--email', required=True, help='Email for job notifications')
    parser.add_argument('--time_per_command', type=int, default=300, help='Estimated time per command in seconds')
    parser.add_argument('--allowed_proteins_file', default=None, help='Path to the file containing the list of allowed proteins')
    parser.add_argument("--rep_lists", default=None, help="Path to the folder containing text files with specific protein lists for each replicate")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    submit_job_script(args)

