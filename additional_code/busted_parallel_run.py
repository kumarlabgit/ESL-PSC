import subprocess
import os
import argparse
from pathlib import Path
import json

def run_command(command):
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        print(output.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error running command {e.cmd}")
        print(e.output.decode())
        exit(1)

def main(alignment_path, tree_path, species_list_path):
    user = os.environ.get('USER')
    if not user:
        print("USER environment variable not set. Exiting.")
        exit(1)
    
    tmp_dir = f"/local_scratch/tmp/{user}/"
    os.makedirs(tmp_dir, exist_ok=True)
    alignment_name = Path(alignment_path).stem
    nex_file_path = os.path.join(tmp_dir, f"{alignment_name}.nex")
    output_path = nex_file_path + ".json"
    
    trim_command = f"/home/user/work/hyphy/HYPHYMP LIBPATH=/home/user/work/hyphy/res /home/user/work/hyphy-analyses/LabelTrees/trim-label-tree.bf --tree {tree_path} --msa {alignment_path} --list {species_list_path} --label FOREGROUND --internal-nodes \"Parsimony\" --output {nex_file_path}"
    print(f"Running trim-label-tree command: {trim_command}")
    run_command(trim_command)
    
    busted_command = f"/home/user/work/hyphy/HYPHYMP LIBPATH=/home/user/work/hyphy/res busted --alignment {nex_file_path} --branches FOREGROUND --output {output_path}" 
    print(f"Running BUSTED command: {busted_command}")
    run_command(busted_command)

    with open(output_path) as json_file:
        data = json.load(json_file)
        p_value = data['test results']['p-value']

    p_value_text_path = f"/local_scratch/tmp/{user}/{alignment_name}_pvalue.txt"

    with open(p_value_text_path, 'w') as pval_file:
        pval_file.write(f"{alignment_name},{p_value}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate running trim-label-tree and BUSTED for a given alignment.")
    parser.add_argument("alignment_path", help="Full path to the alignment file.")
    parser.add_argument("tree_path", help="Full path to the tree file.")
    parser.add_argument("species_list_path", help="Full path to the list of species to label.")
    
    args = parser.parse_args()
    main(args.alignment_path, args.tree_path, args.species_list_path)
