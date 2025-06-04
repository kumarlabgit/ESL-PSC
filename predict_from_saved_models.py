#!/usr/bin/env python3
"""
Apply previously trained ESL-PSC models (new-style .txt files that end in
'_out_feature_weights.txt') to a fresh set of alignments.

Example
-------
python esl_predict_from_models.py \
    --models_dir trained_models/ \
    --prediction_alignments_dir my_alignments/ \
    --response_matrix_path combo_0.txt \
    --output predictions_from_saved_models.csv
"""
import argparse, os, math, re
from collections import defaultdict
from tqdm import tqdm
from Bio import SeqIO
from esl_psc_functions import parse_ESL_weight_line, parse_ESL_weight_label, get_species_to_check

def predict_one_model(model_path, aln_dir, input_species):
    """
    Return (lambda1, lambda2, penalty_term, input_RMSE, species_scores dict)
    for a single saved model file.
    """
    with open(model_path) as fh:
        all_lines = fh.readlines()

    # grab intercept from last line of file
    intercept = parse_ESL_weight_line(all_lines[-1])[1]
    sps = defaultdict(lambda: intercept)

    # infer lambda1, lambda2, penalty term from the file name convention
    m = re.search(r"_l1_([0-9.]+)_l2_([0-9.]+)_out_feature_weights\.txt", model_path)
    lambda1, lambda2 = map(float, m.groups())
    penalty_term = re.search(r'combo_\d+-alignments_(\d+)', model_path)
    penalty_term = penalty_term.group(1) if penalty_term else 'NA'

    for raw in all_lines[:-1]:
        label, weight = parse_ESL_weight_line(raw)
        if weight == 0:
            continue
        gene, pos, aa = parse_ESL_weight_label(label)

        aln_fp = os.path.join(aln_dir, f"{gene}.fas")
        for rec in SeqIO.parse(aln_fp, "fasta"):
            if rec.seq[int(pos)] == aa:
                sps[rec.id] += weight

    # RMSE on input species (if provided)
    rmse = 'NA'
    if input_species:
        err = sum([(sps[sp] - 1) ** 2 if i % 2 == 0 else (sps[sp] + 1) ** 2
                   for i, sp in enumerate(input_species)])
        rmse = math.sqrt(err / len(input_species))

    return lambda1, lambda2, penalty_term, rmse, sps


def main():
    pa = argparse.ArgumentParser(description="Apply saved ESL-PSC models")
    pa.add_argument('--models_dir', required=True,
                    help="Folder with *_out_feature_weights.txt files")
    pa.add_argument('--prediction_alignments_dir', required=True,
                    help="Folder with 2-line fasta alignments for *all* species")
    pa.add_argument('--response_matrix_path',
                    help="Single response matrix (single-matrix runs)")
    pa.add_argument('--response_matrices_dir',
                    help="Directory full of combo_#.txt files (multimatrix runs)")
    pa.add_argument('--output', default="saved_model_predictions.csv",
                    help="Name of the CSV to write")
    args = pa.parse_args()

    # --------- helper for per-model input-species list -----------------------
    def load_input_species(combo_tag: str):
        """
        combo_tag  e.g. 'combo_0'
        returns [] if no matching matrix is available / wanted
        """
        # single-matrix mode
        if args.response_matrix_path and not args.response_matrices_dir:
            return get_species_to_check(args.response_matrix_path, check_order=False)

        # multimatrix mode
        if args.response_matrices_dir:
            path = os.path.join(args.response_matrices_dir, f"{combo_tag}.txt")
            if os.path.exists(path):
                return get_species_to_check(path, check_order=False)
        return []   # no RMSE

    # ------------------------------------------------------------------------
    model_files = sorted(
        f for f in os.listdir(args.models_dir)
        if f.endswith('_out_feature_weights.txt')
    )

    out_lines = ["model_file,lambda1,lambda2,penalty_term,input_RMSE,species,SPS"]

    for mf in tqdm(model_files, smoothing=0, desc="Models"):
        combo_match = re.search(r'(combo_\d+)', mf)
        combo_tag   = combo_match.group(1) if combo_match else None
        input_species = load_input_species(combo_tag) if combo_tag else []

        l1, l2, pen, rmse, sps = predict_one_model(
            os.path.join(args.models_dir, mf),
            args.prediction_alignments_dir,
            input_species
        )

        for sp, score in sps.items():
            if sp in input_species: # skip training species
                continue
            out_lines.append(f"{mf},{l1},{l2},{pen},{rmse},{sp},{score}")

    with open(args.output, 'w') as oh:
        oh.write("\n".join(out_lines))

    print(f"\nPredictions written to: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
