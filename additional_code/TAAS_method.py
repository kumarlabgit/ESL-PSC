import os
import csv
from Bio import SeqIO
import argparse
from tqdm import tqdm

def target_species_specific_substitutions(fasta_dir, species_csv, output_file):
    with open(species_csv, 'r') as f:
        reader = csv.reader(f)
        next(reader)  
        species_data = {row[0]: int(row[1]) for row in reader}

    target_species = {k for k, v in species_data.items() if v == 1}
    non_target_species = {k for k, v in species_data.items() if v == -1}

    substitutions_count = {}

    fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith('.fas')]
    for fasta_file in tqdm(fasta_files, desc="Processing FASTA files"):
        path = os.path.join(fasta_dir, fasta_file)
        sequences = {record.id: str(record.seq) for record in SeqIO.parse(path, 'fasta')}
        
        present_target_species = target_species & set(sequences.keys())
        present_non_target_species = non_target_species & set(sequences.keys())

        if len(present_target_species) < args.min_fg_species:
            print(f"Skipping {fasta_file}: Less than {args.min_fg_species} foreground species present.")
            continue

        substitutions = 0

        for i in range(len(next(iter(sequences.values())))):
            target_aa = {sequences[species][i] for species in present_target_species}
            non_target_aa = {sequences[species][i] for species in present_non_target_species}
                
            if len(target_aa & non_target_aa) == 0:  
                substitutions += 1

        substitutions_count[fasta_file] = substitutions

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Gene', 'Target Species-specific Amino Acid Substitutions'])
        for gene, count in substitutions_count.items():
            writer.writerow([gene, count])

    print(f"Analysis complete! Results saved to {output_file}.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Identify target species-specific amino acid substitutions.")
    parser.add_argument('fasta_dir', type=str, help="Path to directory of fasta files.")
    parser.add_argument('species_csv', type=str, help="Path to csv file with species identifiers.")
    parser.add_argument('--min_fg_species', type=int, default=2, help="Minimum number of foreground species required in an alignment.")

    args = parser.parse_args()

    parent_dir = os.path.dirname(args.fasta_dir)
    output_file = os.path.join(parent_dir, "output.csv")

    target_species_specific_substitutions(args.fasta_dir, args.species_csv, output_file)
