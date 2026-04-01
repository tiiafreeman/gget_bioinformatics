import gget
import pandas as pd
import os
from functools import reduce

matrix_df = pd.read_csv("my_modules/analyzing_motif_matrix/input_folder/SINV_nsP3_motif_matrix.csv", index_col=0)
all_motifs = matrix_df.columns.tolist()

motif_metadata_folder = "my_modules/analyzing_motif_matrix/data/ELM_metadata_folder"
metadata_file_dictionary = {
    "Interactions": {
        "file": os.path.join(motif_metadata_folder, "elm_Interactions.csv"),
        "key": "ELM identifier"
    },
    "GOTerms": {
        "file": os.path.join(motif_metadata_folder, "elm_GOTerms.csv"),
        "key": "ELM"
    },
    "KEGGPathways": {
        "file": os.path.join(motif_metadata_folder, "elm_KEGGPathways.csv"),
        "key": "elm"
    },
    "Classes": {
        "file": os.path.join(motif_metadata_folder, "elm_Classes.csv"),
        "key": "ELMIdentifier"
    },
    "Reactome": {
        "file": os.path.join(motif_metadata_folder, "elm_Reactome.csv"),
        "key": "elm"
    }
}

def metadata_combined_dataframe(metadata_file_dictionary, motif_metadata_folder):
    
    save_path = os.path.join(motif_metadata_folder, "combined_motif_metadata.csv")

    if os.path.exists(save_path):
        print(f"Combined metadata file already exists at {save_path}. Loading existing file.")
        combined_metadata_df = pd.read_csv(save_path)
        
        return combined_metadata_df
    
    else:

        all_metadata_df = []

        for name, info in metadata_file_dictionary.items():
            path = info["file"]

            if os.path.exists(path):
                df = pd.read_csv(path)
                df = df.rename(columns={info["key"]: "Motif"})
                df.columns = [f"{name}_{col}" if col != "Motif" else col for col in df.columns]

                all_metadata_df.append(df)
            else:
                print(f"Warning: File {path} does not exist and will be skipped.")

        combined_metadata_df = reduce(lambda left, right: pd.merge(left, right, on='Motif', how='outer'), all_metadata_df)
        combined_metadata_df = combined_metadata_df.groupby('Motif').agg(
            lambda x: '; '.join(x.dropna().astype(str).unique())
        ).reset_index()

        save_path = os.path.join(motif_metadata_folder, "combined_motif_metadata.csv")
        combined_metadata_df.to_csv(save_path, index=False)

        print(f"Combined metadata saved to {save_path}")

        return combined_metadata_df


combined_metadata_df = metadata_combined_dataframe(metadata_file_dictionary, motif_metadata_folder)


### Crashed from too much data in memory/RAM, make more memory efficient when combining