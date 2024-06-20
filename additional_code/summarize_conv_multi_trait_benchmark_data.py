import argparse
import os
from tqdm import tqdm
import pandas as pd
from pathlib import Path

def read_background(background_dir, combo_name):
    background_csv = os.path.join(background_dir, f"{combo_name}.txt")
    df = pd.read_csv(background_csv)
    df['protein_name'] = df['Gene File'].apply(lambda x: x.split('_')[1])
    df.columns = df.columns.str.replace('"', '')
    return df

def process_ccs(rep_dir, param_dir, background_df):
    ccs_output_file = next((f for f in os.listdir(rep_dir) if f.startswith(f"{param_dir}_") and f.endswith("_CCS_output.csv")), None)
    if ccs_output_file is None:
        return 0, 0

    data = pd.read_csv(os.path.join(rep_dir, ccs_output_file))
    data['protein_name'] = data['Gene File'].apply(lambda x: x.split('_')[1])

    average_ccs_sites = data['Average True Convergence'].mean()

    background_filtered = background_df[~background_df['protein_name'].isin(data['protein_name'])]

    combined = pd.concat([background_filtered, data[['protein_name', 'Average True Convergence']]], ignore_index=True)
    combined.sort_values(by='Average True Convergence', ascending=False, inplace=True)
    top_100_combined = combined.head(100)

    top_100_simulated_count = top_100_combined['protein_name'].isin(data['protein_name']).sum()

    return top_100_simulated_count, average_ccs_sites

def process_esl(rep_dir, param_dir):
    esl_output_file = next((f for f in os.listdir(rep_dir) if f.endswith(f"{param_dir}_gene_ranks.csv")), None)
    if esl_output_file is None:
        return 0

    data = pd.read_csv(os.path.join(rep_dir, esl_output_file))
    top_100_simulated_count = data.head(101)['gene_name'].str.contains('_simulated').sum()

    return top_100_simulated_count

def analyze_method_performance(esl_dir, ccs_dir, background_ccs_dir, species_groups_dir, output_path):
    species_groups_dir_path = Path(species_groups_dir)
    species_combos = [f.stem for f in species_groups_dir_path.iterdir() if f.is_file()]
    results = []

    for combo_name in tqdm(species_combos, smoothing=0):
        background_df = read_background(background_ccs_dir, combo_name)
        ccs_combo_dir = Path(ccs_dir) / combo_name
        esl_combo_dir = Path(esl_dir) / combo_name

        if not ccs_combo_dir.exists() or not esl_combo_dir.exists():
            print(f"Combo directory does not exist for {combo_name}, skipping...")
            continue

        for rep_dir in os.listdir(ccs_combo_dir):
            rep_path = ccs_combo_dir / rep_dir
            if not rep_path.is_dir():
                continue
            for param_dir in os.listdir(rep_path):
                param_path = rep_path / param_dir
                if not param_path.is_dir():
                    continue
                ccs_count, average_ccs_sites = process_ccs(str(param_path), param_dir, background_df)
                esl_count = process_esl(str(esl_combo_dir / rep_dir), param_dir)
                results.append([combo_name, rep_dir, param_dir, ccs_count, average_ccs_sites, esl_count])

    df = pd.DataFrame(results, columns=['Combo', 'Replicate', 'ParamDir', 'CCS_Count', 'Average_CCS_Sites', 'ESL_Count'])
    df.to_csv(output_path, index=False)
    print(f"Analysis summary saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Summarize ESL and CCS results across multiple species sets.")
    parser.add_argument("--esl_dir", required=True, help="Directory for ESL method output data")
    parser.add_argument("--ccs_dir", required=True, help="Directory for CCS method output data")
    parser.add_argument("--background_ccs_dir", required=True, help="Directory with background counts for CCS")
    parser.add_argument("--species_groups_dir", required=True, help="Directory of species groups files for each species set.")
    parser.add_argument("--output_path", required=True, help="Output path for the summary CSV file")
    
    args = parser.parse_args()

    analyze_method_performance(args.esl_dir, args.ccs_dir, args.background_ccs_dir, args.species_groups_dir, args.output_path)

if __name__ == "__main__":
    main()

