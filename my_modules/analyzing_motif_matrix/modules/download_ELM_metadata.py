import pandas as pd
import requests
import io
import os

metadata_folder = "my_modules/analyzing_motif_matrix/data/ELM_metadata_folder"
os.makedirs(metadata_folder, exist_ok=True)

elm_base = "http://elm.eu.org"

elm_endpoints = {
    "Interactions": elm_base + "/infos/browse_elm_interactiondomains.tsv",
    "GOTerms": elm_base + "/goterms.tsv",
    "KEGGPathways": elm_base + "/kegg_terms.tsv",
    "Classes": elm_base + "/elms/elms_index.tsv",
    "Reactome": elm_base + "/reactome_terms.tsv"
}

def download_elm_table(name, url):
    print(f"Downloading {name} metadata...")
    headers = {'User-Agent': 'Mozilla/5.0 (compatible; ELMMetadataDownloader/1.0)'}
    try:
        response = requests.get(url, headers=headers, timeout=60)
        response.raise_for_status()

        raw_text = response.text
        lines = raw_text.splitlines()

        target_headers = ["ELM identifier", "Interaction Domain Id", "Interaction Domain Description", "Interaction Domain Name",
                          "uniprot_id", "elm", "kegg_id", "kegg_name",
                          "uniprot_id","elm","pathway_id", "display_name",
                          "Accession", "GO term", "GO ID"]

        skip_rows = 0
        header_found = False

        for i, line in enumerate(lines):
            if any(h in line for h in target_headers):
                skip_rows = i
                header_found = True
                break

        if not header_found and name == "GOTerms":
            print(f"Warning: Could not find expected headers for {name}. Attempting to read without skipping rows.")
            skip_rows = 5

        df = pd.read_csv(io.StringIO(response.text), sep="\t", skiprows=skip_rows)

        file_path = os.path.join(metadata_folder, f"elm_{name}.csv")
        df.to_csv(file_path, index=False)
        print(f"Saved {name} metadata to {file_path}")
        return df
    except Exception as e:
        print(f"Error downloading {name} metadata: {e}")
        return None

downloaded_data = {}
for name, url in elm_endpoints.items():
    downloaded_data[name] = download_elm_table(name, url)