import pandas as pd
import glob
import os
import argparse

input_folder = "my_modules/ELM_motif_query/output_folder"
output_folder = "my_modules/ELM_motif_query/output_folder"
os.makedirs(output_folder, exist_ok=True)


def create_motif_matrix(input_folder, matrix_name, result_type="labeled"):

    # Select input filename pattern based on result type
    if result_type == "labeled":
        suffix = "_elm_motif_results_labeled.csv"
    elif result_type == "disordered":
        suffix = "_elm_motif_results_disordered.csv"
    elif result_type == "anchor2":
        suffix = "_elm_motif_results_anchor2.csv"

    file_pattern = os.path.join(input_folder, f"*{suffix}")
    files = glob.glob(file_pattern)

    print(f"Looking for: {file_pattern}")
    print(f"Found {len(files)} files.")

    all_data = []

    for path in files:
        file_name = os.path.basename(path).replace(suffix, "")

        df = pd.read_csv(path)

        # Keep one row for each unique ELMIdentifier
        unique_df = df.drop_duplicates(subset=["ELMIdentifier"]).copy()

        print(f"{file_name}: {len(df)} motif hits -> {len(unique_df)} unique ELM identifiers")

        # Save unique version of each input file
        unique_name = f"{file_name}{suffix[:-4]}_unique.csv"
        unique_path = os.path.join(output_folder, unique_name)
        unique_df.to_csv(unique_path, index=False)

        # Add unique motifs to presence/absence matrix
        for motif in unique_df["ELMIdentifier"]:
            all_data.append({
                "Isolate": file_name,
                "Motif": motif,
                "Presence": 1
            })

    if not all_data:
        raise ValueError("No motif data were found. Check the result type and ELMIdentifier column.")

    long_results_df = pd.DataFrame(all_data)

    matrix = long_results_df.pivot(
        index="Isolate",
        columns="Motif",
        values="Presence"
    ).fillna(0).astype(int)

    matrix = matrix.sort_index(axis=1)

    # Save motif matrix
    matrix_file = f"{matrix_name}_{result_type}_motif_matrix.csv"
    matrix.to_csv(os.path.join(output_folder, matrix_file))

    print(f"\nMatrix saved as: {matrix_file}")
    print("\nMotif matrix:")
    print(matrix.head())


if __name__ == "__main__":

    p = argparse.ArgumentParser(description="Create a presence/absence matrix of unique ELM motifs across sequences.")
    p.add_argument("-n", "--name", required=True, help="Name used for the output motif matrix.")
    p.add_argument("-t", "--type", default="labeled", choices=["labeled", "disordered", "anchor2"], help="ELM result type to use. Default: labeled.")
    args = p.parse_args()

    create_motif_matrix(
        input_folder=input_folder,
        matrix_name=args.name,
        result_type=args.type
    )