import os
import csv
import argparse
from Bio import AlignIO
from pathlib import Path
from tqdm import tqdm
import itertools


def get_species_index(alignment, species_id):
    for index, record in enumerate(alignment):
        if record.id == species_id:
            return index
    return None

def count_sites(alignment, outgroup_id, phylogenetic_pairs):
    outgroup_idx = get_species_index(alignment, outgroup_id)
    if outgroup_idx is None:
        return None  

    true_convergence_count = 0
    control_convergence_count = 0

    for conv_id, control_id in phylogenetic_pairs:
        if get_species_index(alignment, conv_id) is None or get_species_index(alignment, control_id) is None:
            return None  

    for i in range(len(alignment[0].seq)):
        outgroup_residue = alignment[outgroup_idx].seq[i]

        if any(alignment[get_species_index(alignment, sp)].seq[i] == '-' for sp in [outgroup_id] + [sp for pair in phylogenetic_pairs for sp in pair]):
            continue

        convergent_residues = [alignment[get_species_index(alignment, conv_id)].seq[i] for conv_id, _ in phylogenetic_pairs]
        control_residues = [alignment[get_species_index(alignment, control_id)].seq[i] for _, control_id in phylogenetic_pairs]

        control_matches_outgroup = all(residue == outgroup_residue for residue in control_residues)

        convergent_residue_counts = {residue: convergent_residues.count(residue) for residue in set(convergent_residues)}
        convergent_residue_counts.pop(outgroup_residue, None)  

        if control_matches_outgroup and max(convergent_residue_counts.values(), default=0) >= 2:
            true_convergence_count += 1

        control_residue_counts = {residue: control_residues.count(residue) for residue in set(control_residues)}
        control_residue_counts.pop(outgroup_residue, None)  

        convergent_matches_outgroup = all(residue == outgroup_residue for residue in convergent_residues)

        if convergent_matches_outgroup and max(control_residue_counts.values(), default=0) >= 2:
            control_convergence_count += 1

    return [len(alignment[0].seq), true_convergence_count, control_convergence_count]

def parse_species_pairs(species_pairs_path):
    with open(species_pairs_path, 'r') as f:
        lines = [line.strip().split(',') for line in f.readlines()]

    return list(itertools.product(*lines))

def analyze_directory(directory, outgroup, species_pairs):
    results = []
    skipped_alignments = 0

    for filename in tqdm(os.listdir(directory), smoothing = 0):
        if filename.endswith(".fas"):
            alignment = AlignIO.read(os.path.join(directory, filename), "fasta")
            site_counts = count_sites(alignment, outgroup, species_pairs)
            if site_counts is None:
                skipped_alignments += 1  
                continue  
            results.append([filename] + site_counts)

    return results, skipped_alignments

def main():
    parser = argparse.ArgumentParser(description='CCS Analysis')
    parser.add_argument('alignment_dir', type=str, help='Path to the directory of alignment files')
    parser.add_argument('outgroup', type=str, help='Species identifier for the outgroup')
    parser.add_argument('species_pairs_file', type=str, help='Path to the txt file with species pairs')
    parser.add_argument('--output_path', type=str, help='path to save output. default is in parent dir of input')
    args = parser.parse_args()

    species_combinations = parse_species_pairs(args.species_pairs_file)

    gene_results = {}

    for combination in species_combinations:
        formatted_combination = list(zip(combination[::2], combination[1::2]))
        results, skipped_alignments = analyze_directory(args.alignment_dir, args.outgroup, formatted_combination)

        for result in results:
            gene_file, num_sites, true_convergence, control_convergence = result
            if gene_file not in gene_results:
                gene_results[gene_file] = {"num_sites": num_sites, "true_totals": [], "control_totals": []}

            gene_results[gene_file]["true_totals"].append(true_convergence)
            gene_results[gene_file]["control_totals"].append(control_convergence)

    directory_name = Path(args.alignment_dir).name
    species_pairs_filename = Path(args.species_pairs_file).stem
    if args.output_path: 
        output_path = args.output_path
    else:
	output_filename = f"{directory_name}_{species_pairs_filename}_CCS_output.csv"
        output_path = Path(args.alignment_dir).parent / output_filename

    csv_data = []
    for gene_file, data in gene_results.items():
        avg_true = sum(data["true_totals"]) / len(data["true_totals"])
        avg_control = sum(data["control_totals"]) / len(data["control_totals"])
        diff = avg_true - avg_control
        csv_data.append([gene_file, data["num_sites"], avg_true, avg_control, diff])

    csv_data.sort(key=lambda x: x[2], reverse=True)

    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Gene File", "Number of Sites", "Average True Convergence", "Average Control Convergence", "Difference True-Control"])
        writer.writerows(csv_data)

    print(f"Number of skipped alignments: {skipped_alignments}")
    print(f"Output file created: {output_path}")

if __name__ == "__main__":
    main()
