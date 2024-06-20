import argparse
import os
import random
from itertools import product

def get_allowed_proteins(file_path):
    if file_path:
        with open(file_path, 'r') as file:
            return set([line.strip() for line in file])
    else:
        return None

def get_replicate_specific_proteins(rep_lists_dir, replicate):
    rep_file_path = os.path.join(rep_lists_dir, f'rep{replicate}_files.txt')
    if os.path.exists(rep_file_path):
        with open(rep_file_path, 'r') as file:
            return [line.strip() for line in file]
    else:
        return None

def main(args):
    os.makedirs(args.output_dir, exist_ok=True)

    allowed_proteins = get_allowed_proteins(args.allowed_proteins_file)

    alignment_files = [f for f in os.listdir(args.alignments_dir) if f.endswith('_NT.fas')]

    if allowed_proteins is not None:
        alignment_files = [f for f in alignment_files if f.split('_')[1] in allowed_proteins]

    with open(os.path.join(args.output_dir, 'simulation_commands.txt'), 'w') as command_file:
        for replicate in range(1, args.num_replicates + 1):
            if args.rep_lists:
                specific_proteins = get_replicate_specific_proteins(args.rep_lists, replicate)
                if specific_proteins:
                    selected_alignments = [f.replace('_simulated', '') for f in specific_proteins]
                else:
                    raise ValueError(f"No specific protein list found for replicate {replicate}. Please check the input.")
            else:
                selected_alignments = random.sample(alignment_files, args.num_alignments_per_replicate)

            for num_sites, scaling_factor in product(args.num_sites_list, args.scaling_factors_list):
                subdir_name = f'rep{replicate}/scale{scaling_factor}_sites{num_sites}'
                subdir_path = os.path.join(args.output_dir, subdir_name)
                os.makedirs(subdir_path, exist_ok=True)

                for alignment in selected_alignments:
                    tree_file = alignment.replace('_NT.fas', '_AA_pruned_tree.nwk')
                    tree_file_path = os.path.join(args.trees_dir, tree_file)
                    alignment_path = os.path.join(args.alignments_dir, alignment)

                    command = (f'python3 csubst_sim_parallel_run.py {alignment_path} {tree_file_path} '
                               f'{subdir_path} {args.foreground_file} --num_sites {num_sites} '
                               f'--foreground_scaling_factor {scaling_factor} '
                               f'--convergent_amino_acids {args.convergent_amino_acids}')
                    command_file.write(command + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate commands for evolutionary simulation runs.")
    parser.add_argument("alignments_dir", help="Path to the folder containing all possible alignments")
    parser.add_argument("num_sites_list", type=lambda s: [int(item) for item in s.split(',')], help="Comma-delimited list of numbers of sites to try (e.g., 3,6,9,12,15)")
    parser.add_argument("scaling_factors_list", type=lambda s: [float(item) for item in s.split(',')], help="Comma-delimited list of foreground scaling factors to try (e.g., 20,40,60,80,100)")
    parser.add_argument("trees_dir", help="Path to the directory containing tree files for simulations")
    parser.add_argument("output_dir", help="Path to the directory to store the output simulations in")
    parser.add_argument("foreground_file", help="Path to the foreground file")
    parser.add_argument("convergent_amino_acids", help="Convergent amino acids parameter for the simulation")
    parser.add_argument("num_alignments_per_replicate", type=int, help="Number of alignments to simulate in each replicate")
    parser.add_argument("num_replicates", type=int, help="Number of replicates to generate")
    parser.add_argument("--allowed_proteins_file", default=None, help="Path to the file containing the list of allowed proteins")
    parser.add_argument("--rep_lists", default=None, help="Path to the folder containing text files with specific protein lists for each replicate")

    args = parser.parse_args()
    main(args)
