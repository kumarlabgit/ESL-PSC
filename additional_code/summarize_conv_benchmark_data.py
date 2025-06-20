import argparse
import os
from tqdm import tqdm
import pandas as pd
import re

def read_background(background_csv, type):
    df = pd.read_csv(background_csv)
    if type == 'csubst':
        df.rename(columns={'Protein': 'protein_name'}, inplace=True)
    elif type == 'ccs':
        df['protein_name'] = df['Gene File'].apply(lambda x: x.split('_')[1])
        df.columns = df.columns.str.replace('"', '')
    print("Columns after processing:", df.columns.tolist())  
    return df

def get_protein_name_from_filename(filename):
    return filename.split('_')[1]

def process_csubst(rep_dir, param_dir, background):
    output_files = [f for f in os.listdir(os.path.join(rep_dir, param_dir)) if f.endswith("_csubst_cb_2.tsv")]
    simulated_data_frames = []
    for f in output_files:
        df = pd.read_csv(os.path.join(rep_dir, param_dir, f), sep='\t')
        protein_name = get_protein_name_from_filename(f)
        df['protein_name'] = protein_name  
        simulated_data_frames.append(df)
    try:
        simulated_data = pd.concat(simulated_data_frames)
    except ValueError as e:
        if str(e) == "No objects to concatenate":
            print(f"Error processing csubst directory {rep_dir}/{param_dir}: No data")
            return 0, 0  
    
    background_filtered = background[~background['protein_name'].isin(simulated_data['protein_name'])]
    combined = pd.concat([background_filtered, simulated_data[['omegaCany2spe', 'OCNCoD', 'protein_name']]], ignore_index=True)
    
    top_100_omega = combined.sort_values(by='omegaCany2spe', ascending=False).head(100)
    top_100_ocn = combined.sort_values(by='OCNCoD', ascending=False).head(100)
    
    count_top_100_omega_simulated = top_100_omega['protein_name'].isin(simulated_data['protein_name']).sum()
    count_top_100_ocn_simulated = top_100_ocn['protein_name'].isin(simulated_data['protein_name']).sum()

    return count_top_100_omega_simulated, count_top_100_ocn_simulated

def process_ccs(rep_dir, param_dir, background):
    ccs_output_file = next((f for f in os.listdir(rep_dir) if f.startswith(param_dir + "_") and f.endswith("_CCS_output.csv")), None)
    if ccs_output_file is None:
        print("no ccs file found for rep, param: ", rep_dir, param_dir)
        return 0,0  

    data = pd.read_csv(os.path.join(rep_dir, ccs_output_file))
    data['protein_name'] = data['Gene File'].apply(lambda x: x.split('_')[1])

    average_ccs_sites = data['Average True Convergence'].mean()

    background_filtered = background[~background['protein_name'].isin(data['protein_name'])]

    data.rename(columns={'Average True Convergence': 'Average True Convergence'}, inplace=True)

    combined = pd.concat([background_filtered, data[['protein_name', 'Average True Convergence']]], ignore_index=True)

    combined = combined.sample(frac=1, random_state=42)  
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

def analyze_method_performance(root_dirs, background_files, output_path):
    csubst_bg = read_background(background_files[0], 'csubst')
    ccs_bg = read_background(background_files[1], 'ccs')
    results = []

    for rep_dir in tqdm([d for d in os.listdir(root_dirs[0]) if os.path.isdir(os.path.join(root_dirs[0], d))], smoothing=0):
        for param_dir in [d for d in os.listdir(os.path.join(root_dirs[0], rep_dir)) if os.path.isdir(os.path.join(root_dirs[0], rep_dir, d))]:
            num_alignments = len([f for f in os.listdir(os.path.join(root_dirs[0], rep_dir, param_dir)) if f.endswith('.fas')])
            
            match = re.match(r'scale(\d+(?:\.\d+)?)_sites(\d+)', param_dir)
            if match:
                scale, sites = match.groups()
                scale = float(scale) if '.' in scale else int(scale)  
                sites = int(sites)  
            else:
                scale, sites = None, None  
            
            csubst_omega, csubst_ocn = process_csubst(root_dirs[0], os.path.join(rep_dir, param_dir), csubst_bg)
            ccs_count, average_ccs_sites = process_ccs(os.path.join(root_dirs[1], rep_dir), param_dir, ccs_bg)
            esl_count = process_esl(os.path.join(root_dirs[2], rep_dir), param_dir)
            
            results.append([rep_dir, scale, sites, num_alignments, csubst_omega, csubst_ocn, ccs_count, esl_count, average_ccs_sites])

    df = pd.DataFrame(results, columns=['Replicate', 'Scale', 'Sites', 'NumAlignments', 'CSUBST_Omega', 'CSUBST_OCN', 'CCS_Count', 'ESL_Count', 'Average_CCS_Sites'])
    df.to_csv(output_path, index=False)
    print(f"Analysis summary saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Analyze performance of convergent evolution detection methods across different parameter combinations.")
    parser.add_argument('csubst_root', help="Root directory for CSUBST method output data")
    parser.add_argument('ccs_root', help="Root directory for CCS method output data")
    parser.add_argument('esl_root', help="Root directory for ESL method output data")
    parser.add_argument('background_csubst', help="CSV file with background values for CSUBST")
    parser.add_argument('background_ccs', help="CSV file with background counts for CCS")
    parser.add_argument('output_path', help="Output path for the summary file")
    
    args = parser.parse_args()

    analyze_method_performance([args.csubst_root, args.ccs_root, args.esl_root], [args.background_csubst, args.background_ccs], args.output_path)

if __name__ == "__main__":
    main()

