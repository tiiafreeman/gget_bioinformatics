import gget
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
import sys
import argparse

sys.path.insert(0, "iupred3")
from iupred3_automation import run_iupred3_on_sequence

input_dir = Path("my_modules/ELM_motif_query/input_folder")
output_dir = Path("my_modules/ELM_motif_query/output_folder")
output_dir.mkdir(parents=True, exist_ok=True)

gget.setup("elm")

def iupred3_filter_motifs_by_disorder(
    sequence_name, sequence,
    regex_df, threshold=0.5,
    iupred_type="long", smoothing="medium", anchor2=False):

    if regex_df is None or regex_df.empty:
        return regex_df, regex_df
    
    out = run_iupred3_on_sequence(
        sequence_name=sequence_name,
        sequence=sequence,
        iupred_type=iupred_type,
        smoothing=smoothing,
        anchor2=anchor2)
    scores = out["scores"]

    # Get amino acid positions of all disordered residues (thresholded as IUPred score > threshold)
    disordered = np.arange(len(scores)+1)[1:][np.array(scores) > threshold]

    # Get amino acid positions of all ordered residues (thresholded as IUPred score < threshold)
    ordered = np.arange(len(scores)+1)[1:][np.array(scores) < threshold]

    def is_subset(arr1, arr2):
        return set(arr2).issubset(set(arr1))

    iupred = []
    for index, row in regex_df.iterrows():
        # Get the start and end amino acid positions of the motif
        start = int(row['motif_start_in_query'])
        end = int(row['motif_end_in_query'])
        motif_positions = range(start, end + 1)

        # Check if all positions covered by the motif are within ordered or disordered IUPred scores thresholds
        if is_subset(ordered, list(motif_positions)):
            iupred.append("ordered")

        elif is_subset(disordered, list(motif_positions)):
            iupred.append("disordered")

        else:
            iupred.append("inbetween")

    labeled_df = regex_df.copy()
    labeled_df["IUPred"] = iupred

    disordered_only_df = labeled_df[labeled_df["IUPred"] == "disordered"].copy()

    if anchor2:
        anchor2_scores = out["anchor2"]

        anchor2_positive = np.arange(len(anchor2_scores) + 1)[1:][
            np.array(anchor2_scores) > threshold
        ]

        anchor2_labels = []

        for index, row in disordered_only_df.iterrows():
            start = int(row["motif_start_in_query"])
            end = int(row["motif_end_in_query"])
            motif_positions = range(start, end + 1)

            if is_subset(anchor2_positive, list(motif_positions)):
                anchor2_labels.append("positive")
            else:
                anchor2_labels.append("negative")

        disordered_only_df["ANCHOR2"] = anchor2_labels

        anchor2_only_df = disordered_only_df[disordered_only_df["ANCHOR2"] == "positive"].copy()

    labeled_name = output_dir / f"{sequence_name}_elm_motif_results_labeled.csv"
    disordered_name = output_dir / f"{sequence_name}_elm_motif_results_disordered.csv"

    labeled_df.to_csv(labeled_name, index=False)
    disordered_only_df.to_csv(disordered_name, index=False)

    if anchor2:
        anchor2_name = output_dir / f"{sequence_name}_elm_motif_results_anchor2.csv"
        anchor2_only_df.to_csv(anchor2_name, index=False)

    return labeled_df, disordered_only_df, anchor2_only_df


def search_elm_motifs(input_dir):

    pd.set_option('display.max_columns', None)

    for input_file in input_dir.glob("*.fasta"):
        for record in SeqIO.parse(input_file, "fasta"):
            print(f"Processing {record.id}...")

            ortholog_df, regex_df = gget.elm(str(record.seq), uniprot=False, expand=True)

            if ortholog_df is not None and not ortholog_df.empty:
                ortholog_name = output_dir / f"{record.id}_elm_ortholog_results.csv"

                # ortholog df contains info about DIAMOND alignment
                ortholog_df.to_csv(ortholog_name, index=False)

            if regex_df is not None and not regex_df.empty:

                if args.filter:
                    labeled_df, disordered_only_df, anchor2_only_df = iupred3_filter_motifs_by_disorder(
                        sequence_name=record.id,
                        sequence=str(record.seq),
                        regex_df=regex_df,
                        threshold=0.5,
                        iupred_type=args.iupred_type,
                        smoothing=args.smoothing,
                        anchor2=args.anchor2,
                    )
                else:
                    regex_name = output_dir / f"{record.id}_elm_motif_results.csv"
                    regex_df.to_csv(regex_name, index=False)

    return ortholog_df, regex_df

if __name__ == "__main__":

    p = argparse.ArgumentParser(description="Search for ELM motifs in FASTA files and filter by IUPred3 disorder scores.")
    p.add_argument("--filter", action="store_true", help="Filter ELM motifs to disordered regions using IUPred3.")
    p.add_argument("-i", "--iupred_type", default="long", choices=["long", "short"], help="Type of IUPred3 prediction to run (long or short).")
    p.add_argument("-s", "--smoothing", default="medium", choices=["none", "medium", "strong"], help="Smoothing level for IUPred3 scores.")
    p.add_argument("-a", "--anchor2", default=False, action="store_true", help="Use anchor2 filtering.")
    args = p.parse_args()


    search_elm_motifs(input_dir)