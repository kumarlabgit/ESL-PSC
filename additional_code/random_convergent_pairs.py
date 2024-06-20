import os
import argparse
from ete3 import Tree
import random
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description="Generate random species sets for convergent evolution studies.")
    parser.add_argument("--tree_file", type=str, required=True, help="Path to the Newick tree file.")
    parser.add_argument("--num_combinations", type=int, required=True, help="Number of combinations to select.")
    parser.add_argument("--num_species_in_combination", type=int, required=True, help="Number of species in each combination.")
    parser.add_argument("--exclude_species_file", type=str, required=True, help="Path to the file containing species to exclude.")
    parser.add_argument("--output_dir_pairs", type=str, required=True, help="Path to the output directory for species pairs.")
    parser.add_argument("--output_dir_convergent", type=str, required=True, help="Path to the output directory for convergent species.")
    parser.add_argument("--min_distance", type=float, required=True, help="Minimum evolutionary distance between the convergent and control species.")
    return parser.parse_args()

def read_exclude_species(file_path):
    """Read species to exclude from a file."""
    try:
        with open(file_path, 'r') as file:
            return [line.strip() for line in file if line.strip()]
    except FileNotFoundError:
        raise ValueError(f"Excluded species file not found: {file_path}")
    except Exception as e:
        raise ValueError(f"Error reading excluded species file: {e}")

def load_tree(tree_file):
    try:
        return Tree(tree_file)
    except Exception as e:
        raise ValueError(f"Failed to load the tree from {tree_file}: {e}")

def get_valid_species(tree, exclude_species):
    """Return a list of species excluding the specified ones and their descendants."""
    valid_species = [leaf.name for leaf in tree.iter_leaves() if leaf.name not in exclude_species]
    return set(valid_species)

def find_valid_control(convergent_species_node, tree, exclude_species_nodes, min_distance):
    """
    Find a valid control species by collecting up to three candidates
    that respect the exclusion criteria, ensuring control is not the same
    as the convergent species, and then randomly selecting one.
    """
    potential_controls = []
    species_node = convergent_species_node  

    while species_node.up:
        parent = species_node.up
        siblings = [s for s in parent.get_descendants() if s not in exclude_species_nodes and s.is_leaf() and s.name.strip()]
        random.shuffle(siblings)  

        for sibling in siblings:
            if sibling != convergent_species_node and convergent_species_node.get_distance(sibling) >= min_distance:
                mrca_node = tree.get_common_ancestor(convergent_species_node, sibling)
                if not any(descendant in exclude_species_nodes for descendant in mrca_node.iter_descendants()):
                    potential_controls.append(sibling)
                    if len(potential_controls) == 3:
                        break
        if len(potential_controls) == 3:
            break  

        species_node = parent  

    return random.choice(potential_controls) if potential_controls else None

def update_exclusions(mrca_node, exclude_species_nodes):
    """
    Update the exclusion set to include the MRCA and its descendants.
    """
    descendants = {descendant for descendant in mrca_node.iter_descendants()}
    exclude_species_nodes.update(descendants)
    exclude_species_nodes.add(mrca_node)

def select_species_pairs(tree, num_species, exclude_species_nodes, min_distance):
    """
    Select species pairs ensuring that exclusions are respected.
    """
    def attempt_selection():
        valid_species_nodes = {leaf for leaf in tree.iter_leaves() if leaf not in exclude_species_nodes}
        selected_pairs = []
        
        while len(selected_pairs) * 2 < num_species and valid_species_nodes:
            convergent_species_node = random.choice(list(valid_species_nodes))
            valid_control_node = find_valid_control(convergent_species_node, tree, exclude_species_nodes, min_distance)
            if valid_control_node:
                selected_pairs.append((convergent_species_node, valid_control_node))
                mrca_node = tree.get_common_ancestor(convergent_species_node, valid_control_node)
                update_exclusions(mrca_node, exclude_species_nodes)  
                valid_species_nodes.difference_update(exclude_species_nodes)  
            else:
                valid_species_nodes.remove(convergent_species_node)
        
        return selected_pairs if len(selected_pairs) * 2 >= num_species else None

    trials = 0
    result = None
    while result is None and trials < 100:
        result = attempt_selection()
        trials += 1
        if trials >= 100:
            raise ValueError("Failed to find enough valid species pairs after 100 attempts. Please check your input criteria and try again.")
    return result

def write_output_files(selected_pairs, output_dir_pairs, output_dir_convergent, combination_index):
    os.makedirs(output_dir_pairs, exist_ok=True)
    os.makedirs(output_dir_convergent, exist_ok=True)

    pairs_file_path = os.path.join(output_dir_pairs, f"combination_{combination_index}.txt")
    convergent_file_path = os.path.join(output_dir_convergent, f"combination_{combination_index}.txt")

    with open(pairs_file_path, 'w') as f_pairs:
        for species_node, control_node in selected_pairs:
            species_name = species_node.name.replace('--', '')
            control_name = control_node.name.replace('--', '')
            f_pairs.write(f"{species_name}\n{control_name}\n")

    with open(convergent_file_path, 'w') as f_convergent:
        for i, (species_node, _) in enumerate(selected_pairs, start=1):
            species_name = species_node.name.replace('--', '').strip()
            f_convergent.write(f"{i}\t{species_name}\n")

def write_output_files(selected_pairs, output_dir_pairs, output_dir_convergent, combination_index):
    os.makedirs(output_dir_pairs, exist_ok=True)
    os.makedirs(output_dir_convergent, exist_ok=True)

    pairs_file_path = os.path.join(output_dir_pairs, f"combination_{combination_index}.txt")
    convergent_file_path = os.path.join(output_dir_convergent, f"combination_{combination_index}.txt")

    with open(pairs_file_path, 'w') as f_pairs:
        for species_node, control_node in selected_pairs:
            species_name = species_node.name.replace('--', '')
            control_name = control_node.name.replace('--', '')
            f_pairs.write(f"{species_name}\n{control_name}\n")

    with open(convergent_file_path, 'w') as f_convergent:
        for i, (species_node, _) in enumerate(selected_pairs, start=1):
            species_name = species_node.name.replace('--', '').strip()
            f_convergent.write(f"{i}\t{species_name}\n")

def main():
    args = parse_args()
    tree = load_tree(args.tree_file)

    os.makedirs(args.output_dir_pairs, exist_ok=True)
    os.makedirs(args.output_dir_convergent, exist_ok=True)

    excluded_species_names = read_exclude_species(args.exclude_species_file)
    max_retries = 100  

    for i in tqdm(range(args.num_combinations), smoothing=0):
        success = False
        for retry in range(max_retries):
            try:
                excluded_species_nodes = {tree & name for name in excluded_species_names if name in tree}
                selected_pairs = select_species_pairs(tree, args.num_species_in_combination, excluded_species_nodes, args.min_distance)
                write_output_files(selected_pairs, args.output_dir_pairs, args.output_dir_convergent, i + 1)
                success = True
                break  
            except ValueError as e:
                print(f"Attempt {retry + 1}: {e}")

        if not success:
            raise RuntimeError(f"Failed to generate valid species pairs for combination {i + 1} after {max_retries} retries.")

    print(f"Output files generated in {args.output_dir_pairs} and {args.output_dir_convergent}")

if __name__ == "__main__":
    main()

