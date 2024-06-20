import argparse
import subprocess
import os
from pathlib import Path

def submit_individual_job_scripts(args):
    allowed_proteins_arg = f"--allowed_proteins_file {args.allowed_proteins_file}" if args.allowed_proteins_file else ""
    rep_lists_arg = f"--rep_lists {args.rep_lists}" if args.rep_lists else ""

    if args.foreground_dir:
        for fg_file in os.listdir(args.foreground_dir):
            fg_path = os.path.join(args.foreground_dir, fg_file)
            sub_output_dir = os.path.join(args.output_dir, Path(fg_file).stem)
            os.makedirs(sub_output_dir, exist_ok=True)
            generate_and_submit_job_script(fg_path, sub_output_dir, args, allowed_proteins_arg, rep_lists_arg)
    else:
        generate_and_submit_job_script(args.foreground_file, args.output_dir, args, allowed_proteins_arg, rep_lists_arg)

def generate_and_submit_job_script(fg_path, output_dir, args, allowed_proteins_arg, rep_lists_arg):
    combo_name = Path(fg_path).stem
    job_name = f"{combo_name}_{os.path.basename(args.output_dir.rstrip('/'))}"
    script_filename = f"{job_name}_job_script.sh"
    
    job_script_content = f"""#!/bin/sh
#PBS -l walltime={args.walltime}
#PBS -N {job_name}
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -m bae
#PBS -M {args.email}


cd $PBS_O_WORKDIR
module load python/3.8.5
source ~/env_name/bin/activate

python3 sim_preprocess_for_esl_benchmarks.py {args.alignments_dir} {','.join(map(str, args.num_sites))} {','.join(map(str, args.scaling_factors))} {args.trees_dir} {output_dir} {fg_path} {args.convergent_aa} {args.num_align_per_rep} {args.num_replicates} {allowed_proteins_arg} {rep_lists_arg}

torque-launch -p {combo_name}_chkpt.txt {os.path.join(output_dir, 'simulation_commands.txt')}
"""

    with open(script_filename, "w") as file:
        file.write(job_script_content)

    subprocess.run(["qsub", script_filename])

def parse_arguments():
    parser = argparse.ArgumentParser(description='HPC Job Submission Script for running csubst simulations of convergent evolution')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--foreground_file', help='Foreground file')
    group.add_argument('--foreground_dir', help='Directory containing all foreground files')
    parser.add_argument('--alignments_dir', required=True, help='Path to alignments directory')
    parser.add_argument('--num_sites', type=lambda s: [int(item) for item in s.split(',')], required=True, help='List of numbers of sites (comma-separated)')
    parser.add_argument('--scaling_factors', type=lambda s: [float(item) for item in s.split(',')], required=True, help='List of scaling factors (comma-separated)')
    parser.add_argument('--trees_dir', required=True, help='Path to trees directory')
    parser.add_argument('--output_dir', required=True, help='Path to output directory')
    parser.add_argument('--convergent_aa', required=True, help='Convergent amino acids')
    parser.add_argument('--num_align_per_rep', type=int, required=True, help='Number of alignments per replicate')
    parser.add_argument('--num_replicates', type=int, required=True, help='Number of replicates')
    parser.add_argument('--email', required=True, help='Email for job notifications')
    parser.add_argument('--walltime', required=True, help='Walltime for the job')
    parser.add_argument('--nodes', required=True, help='Number of nodes for the job')
    parser.add_argument('--allowed_proteins_file', default=None, help='Path to the file containing the list of allowed proteins')
    parser.add_argument("--rep_lists", default=None, help="Path to the folder containing text files with specific protein lists for each replicate, for reproducing")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    submit_individual_job_scripts(args)




