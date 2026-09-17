import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import argparse

input_folder = "my_modules/ELM_motif_query/output_folder"
output_folder = "my_modules/ELM_motif_query/output_folder"
os.makedirs(output_folder, exist_ok=True)


def motif_conservation(input_file=None):

    # Use specified matrix or find motif matrix in default folder
    if input_file:
        matrix_file = input_file

    else:
        matrix_files = glob.glob(os.path.join(input_folder, "*_motif_matrix.csv"))

        if len(matrix_files) == 0:
            raise FileNotFoundError(f"No file ending in '_motif_matrix.csv' found in {input_folder}")

        if len(matrix_files) > 1:
            raise ValueError(f"Multiple motif matrix files found in {input_folder}: {matrix_files}")

        matrix_file = matrix_files[0]

    print(f"Using motif matrix: {matrix_file}")

    matrix_df = pd.read_csv(matrix_file, index_col=0)

    # Calculate percentage of sequences containing each motif
    motif_percentages = (matrix_df.sum() / len(matrix_df)) * 100

    stats_df = pd.DataFrame({
        "Motif": motif_percentages.index,
        "Percentage": motif_percentages.values
    })

    # Group motifs by ELM class prefix
    stats_df["Prefix"] = stats_df["Motif"].str.split("_").str[0]

    for prefix, group in stats_df.groupby("Prefix"):

        plt.figure(figsize=(10, 6))

        group = group.sort_values("Percentage", ascending=False)

        plt.bar(group["Motif"], group["Percentage"])

        plt.title(f"Percentage of Sequences with {prefix} Motifs")
        plt.xlabel("Motif ID")
        plt.ylabel("Percentage of Sequences")
        plt.xticks(rotation=45, ha="right")
        plt.ylim(0, 105)

        for i, v in enumerate(group["Percentage"]):
            plt.text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=9)

        plt.tight_layout()

        save_path = os.path.join(
            output_folder,
            f"{prefix}_motif_conservation.png"
        )

        plt.savefig(save_path)
        print(f"Saved plot for {prefix} motifs to {save_path}")
        plt.close()


if __name__ == "__main__":

    p = argparse.ArgumentParser(description="Plot the percentage of viral sequences containing each ELM motif.")
    p.add_argument("-i", "--input", default=None, help="Path to a motif matrix CSV. Default: automatically use *_motif_matrix.csv from the ELM output folder.")
    args = p.parse_args()

    motif_conservation(input_file=args.input)