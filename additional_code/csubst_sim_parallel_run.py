import argparse
import os
import subprocess
import tempfile
from Bio import SeqIO, Seq

def run_csubst_simulation(alignment_path, tree_file_path, output_folder, foreground_file, args):
    temp_dir_base = os.getenv('TMPDIR', tempfile.gettempdir())
    
    with tempfile.TemporaryDirectory(dir=temp_dir_base) as temp_dir:
        command = [
            "csubst", "simulate",
            "--iqtree_exe", "iqtree2",
            "--alignment_file", alignment_path,
            "--rooted_tree_file", tree_file_path,
            "--iqtree_redo", "no",
            "--num_simulated_site", str(args.num_sites),
            "--percent_convergent_site", "100",  
            "--foreground", foreground_file,
            "--optimized_branch_length", "yes",
            "--foreground_scaling_factor", str(args.foreground_scaling_factor),
            "--foreground_omega", str(args.foreground_omega),
            "--convergent_amino_acids", args.convergent_amino_acids
        ]

        subprocess.run(command, check=True, cwd=temp_dir)

        simulated_alignment_path = os.path.join(temp_dir, "simulate.fa")
        simulated_alignment_dict = SeqIO.to_dict(SeqIO.parse(simulated_alignment_path, "fasta"))
        
        original_alignment_dict = SeqIO.to_dict(SeqIO.parse(alignment_path, "fasta"))
        for seq_id, original_record in original_alignment_dict.items():
            simulated_record = simulated_alignment_dict.get(seq_id)
            if simulated_record:
                original_seq_str = str(original_record.seq)
                simulated_seq_str = str(simulated_record.seq)
                combined_seq = simulated_seq_str + original_seq_str[args.num_sites * 3:]
                original_record.seq = Seq.Seq(combined_seq)

        combined_alignment_path = os.path.join(output_folder, os.path.basename(alignment_path).replace('.fas', '_simulated.fas'))
        SeqIO.write(original_alignment_dict.values(), combined_alignment_path, "fasta-2line")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a CSUBST evolutionary simulation on a single alignment.")
    parser.add_argument("alignment_file", help="Path to the alignment file")
    parser.add_argument("tree_file", help="Path to the tree file")
    parser.add_argument("output_folder", help="Path to the output folder")
    parser.add_argument("foreground_file", help="Path to the foreground file")
    parser.add_argument("--num_sites", type=int, required=True, help="Exact number of codons to replace in the original alignment")
    parser.add_argument("--foreground_scaling_factor", type=float, default=2, help="Foreground scaling factor")
    parser.add_argument("--foreground_omega", type=float, default=5, help="Foreground omega")
    parser.add_argument("--convergent_amino_acids", default="random1", help="Convergent amino acids")

    args = parser.parse_args()
    run_csubst_simulation(args.alignment_file, args.tree_file, args.output_folder, args.foreground_file, args)

