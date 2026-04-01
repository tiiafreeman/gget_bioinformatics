import pandas as pd
import glob
import os

input_folder = "my_modules/ELM_motif_query/input_folder"
output_folder = "my_modules/ELM_motif_query/output_folder"
os.makedirs(output_folder, exist_ok=True)

file_pattern = input_folder + "/*_translation_elm_motif_results_labeled_unique.csv"

all_data = []

for path in glob.glob(file_pattern):
    file_name = os.path.basename(path).replace("_translation_elm_motif_results_labeled_unique.csv", "")
    df = pd.read_csv(path)

    motifs = df['ELMIdentifier'].unique()

    for motif in motifs:
        all_data.append({
            'Isolate': file_name,
            'Motif': motif,
            'Presence': 1
        })

long_results_df = pd.DataFrame(all_data)
matrix = long_results_df.pivot(index='Isolate', columns='Motif', values='Presence').fillna(0).astype(int)
matrix = matrix.sort_index(axis=1)

matrix.to_csv(os.path.join(output_folder, "SINV_nsP3_motif_matrix.csv"))
print(matrix.head())