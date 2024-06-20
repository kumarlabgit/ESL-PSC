import argparse
import os
from Bio import Phylo, SeqIO
from tqdm import tqdm
import csv

def get_common_ancestor(tree, species):
    return tree.common_ancestor(species)

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

def get_sequence_by_id(file_path, identifier, log_file_path):
    """Fetch a sequence from a fasta file based on an identifier."""
    for record in SeqIO.parse(file_path, 'fasta'):
        if record.id == identifier:
            return record.seq
    write_log(f"Warning: Sequence not found for ID {identifier} in {file_path}", log_file_path)
    return None

def get_clade_from_name(tree, name):
    """Get Clade object from tree based on species name."""
    for terminal in tree.get_terminals():
        if terminal.name == name:
            return terminal
    return None

def write_log(message, log_file_path):
    """Write a message to the log file."""
    with open(log_file_path, 'a') as log_file:  
        log_file.write(message + '\n')

def main():
    parser = argparse.ArgumentParser(description="Find convergent amino acid residues.")
    parser.add_argument("protein_dir", help="Path to protein sequence alignments")
    parser.add_argument("ancestral_dir", help="Path to ancestral sequence alignments")
    parser.add_argument("species_file", help="TXT file with comma-delimited lists of species identifiers on each line. pairs of lines (lists) denote separate comparisons to run")
    parser.add_argument("--protein_suffix", default="_AA.fas", help="Suffix for protein alignment files")
    parser.add_argument("--ancestral_suffix", default="_iqtree_ancestral_sequences.fasta", help="Suffix for ancestral alignments")
    parser.add_argument("--tree_suffix", default="_slac_tree.nwk", help="Suffix for tree files")
    args = parser.parse_args()

    with open(args.species_file, 'r') as f:
        lines = [line.strip().split(',') for line in f]
        species_pairs = [(lines[i], lines[i+1]) for i in range(0, len(lines), 2)]

    results = []

    log_file_path = os.path.join(os.path.dirname(args.ancestral_dir), "convergence_checker_log_file.txt")

    ancestral_dir_name = os.path.basename(os.path.normpath(args.ancestral_dir))
    species_list_name = os.path.splitext(os.path.basename(args.species_file))[0]
    output_file_name = f"{ancestral_dir_name}_{species_list_name}_convergence_results.csv"
    output_file_path = os.path.join(os.path.dirname(args.ancestral_dir), output_file_name)

    with open(log_file_path, "w") as _:
        pass
    
    for ancestral_file in tqdm(os.listdir(args.ancestral_dir)):
        if not ancestral_file.endswith(args.ancestral_suffix):
            continue
        basename = ancestral_file[:-len(args.ancestral_suffix)]

        protein_file = os.path.join(args.protein_dir, basename + args.protein_suffix)
        tree_file = os.path.join(args.ancestral_dir, basename + args.tree_suffix)

        if not os.path.exists(protein_file) or not os.path.exists(tree_file):
            with open(log_file_path, "a") as log:
                log.write(f"Missing file for {basename}\n")
            continue

        tree = Phylo.read(tree_file, 'newick')

        convergence_counts = []
        for species_list_a, species_list_b in species_pairs:
            terminal_species = [clade.name for clade in tree.get_terminals()]

            species_list_a = [species for species in species_list_a if species in terminal_species]
            species_list_b = [species for species in species_list_b if species in terminal_species]

            if not species_list_a or not species_list_b:
                with open(log_file_path, "a") as log:
                    log.write(f"Either {species_list_a} or {species_list_b} do not have any species in tree for {basename}\n")
                convergence_counts.append('')
                continue

            if len(species_list_a) == 1:
                common_ancestor_a_name = species_list_a[0]
                common_ancestor_a = get_clade_from_name(tree, common_ancestor_a_name)
            else:
                common_ancestor_a = get_common_ancestor(tree, species_list_a)
                
            if len(species_list_b) == 1:
                common_ancestor_b_name = species_list_b[0]
                common_ancestor_b = get_clade_from_name(tree, common_ancestor_b_name)
            else:
                common_ancestor_b = get_common_ancestor(tree, species_list_b)

            parent_a = get_parent(tree, common_ancestor_a)
            parent_b = get_parent(tree, common_ancestor_b)

            seq_a = get_sequence_by_id(os.path.join(args.ancestral_dir, ancestral_file), common_ancestor_a.name, log_file_path) if len(species_list_a) > 1 else get_sequence_by_id(protein_file, species_list_a[0], log_file_path)
            seq_b = get_sequence_by_id(os.path.join(args.ancestral_dir, ancestral_file), common_ancestor_b.name, log_file_path) if len(species_list_b) > 1 else get_sequence_by_id(protein_file, species_list_b[0], log_file_path)
            seq_parent_a = get_sequence_by_id(os.path.join(args.ancestral_dir, ancestral_file), parent_a.name, log_file_path)
            seq_parent_b = get_sequence_by_id(os.path.join(args.ancestral_dir, ancestral_file), parent_b.name, log_file_path)

            if None in [seq_a, seq_b, seq_parent_a, seq_parent_b]:
                write_log(f"Missing ancestral or protein sequence for {ancestral_file} with species pairs {species_list_a} and {species_list_b}", log_file_path)
                continue
            
            count = 0
            write_log(f"Sequences for {basename}:", log_file_path)
            write_log(f"seq_a ({common_ancestor_a.name}): {seq_a}", log_file_path)
            write_log(f"seq_b ({common_ancestor_b.name}): {seq_b}", log_file_path)
            write_log(f"seq_parent_a ({parent_a.name}): {seq_parent_a}", log_file_path)
            write_log(f"seq_parent_b ({parent_b.name}): {seq_parent_b}", log_file_path)
            for (res_a, res_b, res_pa, res_pb) in zip(seq_a, seq_b, seq_parent_a, seq_parent_b):
                if '-' in (res_a, res_b, res_pa, res_pb):
                    continue  

                if res_a == res_b and res_a != res_pa and res_b != res_pb:
                    count += 1

            convergence_counts.append(count)

        protein_length = len(next(SeqIO.parse(protein_file, 'fasta')).seq)
        results.append([basename, protein_length] + convergence_counts)

    with open(output_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Protein', 'Length'] + ['&'.join(pair[0]) + '/' + '&'.join(pair[1]) for pair in species_pairs])
        writer.writerows(results)

    print(f"Output saved to: {output_file_path}")
    print(f"Log saved to: {log_file_path}")

if __name__ == "__main__":
    main()

