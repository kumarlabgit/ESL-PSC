import subprocess
from tqdm import tqdm
import argparse
import os
import tempfile

def run_csubst_simulations(args):
    output_folder_name = "simulated_alignments" + args.output_suffix
    output_folder_path = os.path.join(os.getcwd(), output_folder_name)
    os.makedirs(output_folder_path, exist_ok=True)

    log_file_path = os.path.join(output_folder_path, "csubst_errors.log")

    command_template = [
        "csubst", "simulate",
        "--iqtree_exe", "iqtree2",
        "--alignment_file", args.alignment_file,
        "--rooted_tree_file", args.tree_file,
        "--iqtree_redo", "no",
        "--foreground", args.foreground_file,
        "--optimized_branch_length", "yes",
        "--percent_biased_sub", "99",
        "--num_simulated_site", str(args.num_simulated_site)
    ]

    convergent_params = [
        "--foreground_scaling_factor", str(args.foreground_scaling_factor),
        "--foreground_omega", str(args.foreground_omega),
        "--convergent_amino_acids", args.convergent_amino_acids,
        "--percent_convergent_site", str(args.percent_convergent_site),
    ]
    neutral_params = [
        "--foreground_scaling_factor", "1",
        "--foreground_omega", ".2",
        "--convergent_amino_acids", "random0",
    ]

    print("running simulations")

    with tqdm(total=args.num_conv_alignments + args.num_neutral_alignments, desc="Running simulations") as progress:
        for i in range(args.num_conv_alignments):
            unique_filename = f"convergent_{i + 1}.fas"
            run_simulation(command_template + convergent_params, log_file_path, output_folder_path, unique_filename)
            progress.update(1)
        for i in range(args.num_neutral_alignments):
            unique_filename = f"neutral_{i + 1}.fas"
            run_simulation(command_template + neutral_params, log_file_path, output_folder_path, unique_filename)
            progress.update(1)
    
    print("Output folder: " + output_folder_path)

def run_simulation(command, log_file_path, output_folder_path, unique_filename):
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=temp_dir)
            simulated_file_path = os.path.join(temp_dir, "simulate.fa")  
            unique_file_path = os.path.join(output_folder_path, unique_filename)
            os.rename(simulated_file_path, unique_file_path)
        except subprocess.CalledProcessError:
            with open(log_file_path, 'a') as log_file:
                log_file.write(f"Error with simulation: {' '.join(command)}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CSUBST evolutionary simulations on a single alignment.")
    parser.add_argument("alignment_file", help="Path to the alignment file")
    parser.add_argument("tree_file", help="Path to the tree file")
    parser.add_argument("foreground_file", help="Path to the foreground file")
    parser.add_argument("--num_simulated_site", type=int, default=500, help="Number of sites to be simulated (default 500)")
    parser.add_argument("--num_conv_alignments", type=int, default=10, help="Number of alignments with convergent sites (default 10)")
    parser.add_argument("--num_neutral_alignments", type=int, default=90, help="Number of neutral alignments (default 90)")
    parser.add_argument("--foreground_scaling_factor", type=float, default=2, help="Foreground scaling factor for convergent genes")
    parser.add_argument("--foreground_omega", type=float, default=5, help="Foreground omega for convergent genes")
    parser.add_argument("--convergent_amino_acids", default="random1", help="Convergent amino acids for convergent genes")
    parser.add_argument("--percent_convergent_site", type=int, default=5, help="Percent convergent site for convergent genes")
    parser.add_argument("--output_suffix", default="", help="Optional suffix to append to the output folder name")
    
    args = parser.parse_args()
    run_csubst_simulations(args)
