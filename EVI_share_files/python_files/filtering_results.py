from pathlib import Path
import pandas as pd

input_dir = Path("my_modules/ELM_motif_query/input_folder")
output_dir = Path("my_modules/ELM_motif_query/output_folder")
output_dir.mkdir(parents=True, exist_ok=True)

columns_to_keep = [
    "ELMIdentifier",
    "FunctionalSiteName",
    "ELMType",
    "Description",
    "InteractionDomainDescription",
    "InteractionDomainName",
    "Regex",
    "Instances (Matched Sequence)",
    "motif_start_in_query",
    "motif_end_in_query"
]

for file in input_dir.glob("*labeled.csv"):
    print(f"Processing {file}...")

    df = pd.read_csv(file)

    missing_columns = [c for c in columns_to_keep if c not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in input CSV: {missing_columns}")

    df_unique = (
        df[columns_to_keep]
        .drop_duplicates(subset=columns_to_keep, keep="first")
    )

    output_csv = output_dir / f"{file.stem}_unique.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df_unique.to_csv(output_csv, index=False)