import pandas as pd
from tqdm import tqdm
import requests
import time
from chembl_webresource_client.new_client import new_client


df = pd.read_csv("inbetween_data/matched_reactions_deduplicated.csv")
df = df.drop_duplicates(subset=["uniprot_id"])

print("Unique UniProt IDs:", len(df))

uniprot_id_string = df["uniprot_id"].str.cat(sep=",")

# Map UniProt IDs to ChEMBL target IDs using UniProt's ID mapping API
url = "https://rest.uniprot.org/idmapping/run"
data = {"ids": uniprot_id_string, "from": "UniProtKB_AC-ID", "to": "ChEMBL"}

response = requests.post(url, data=data)
response.raise_for_status()
job_id_response = response.json()
job_id = job_id_response.get("jobId")

if job_id is not None:
    print("Job ID:", job_id)
else:
    print("Failed to retrieve Job ID")
    exit(1)

time.sleep(5)

url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
response = requests.get(url)
response.raise_for_status()
status_data = response.json()
print("Job status:", status_data.get("jobStatus"))

url = f"https://rest.uniprot.org/idmapping/stream/{job_id}?format=json"
response = requests.get(url, stream=True)
response.raise_for_status()
results_data = response.json()
mappings = results_data.get("results", [])
print(f"Successfully retrieved {len(mappings)} mappings.")

mapping_df = pd.DataFrame(mappings)
mapping_df = mapping_df.rename(columns={"from": "uniprot_id", "to": "chembl_id"})
mapping_df = mapping_df[["uniprot_id", "chembl_id"]]
mapping_df = pd.merge(df, mapping_df, on="uniprot_id", how="left")
mapping_df = mapping_df[["uniprot_id", "chembl_id"]]

mapping_df_no_duplicates = mapping_df.drop_duplicates(subset=["uniprot_id"])
if len(mapping_df_no_duplicates) != len(mapping_df):
    print("Warning: Duplicate UniProt IDs were dropped. Look them up!")
    dropped_indices = mapping_df.index.difference(mapping_df_no_duplicates.index)
    print(f"Number of dropped: {len(dropped_indices)}")
    dropped_rows = mapping_df.loc[dropped_indices]
    print("Dropped rows:")
    print(dropped_rows)

mapping_df_no_duplicates.to_csv("results_target_class.csv", index=False)


df = pd.read_csv("results_target_class.csv")

# Store target classification descriptions
target_classification = []
for index, row in tqdm(
    df.iterrows(),
    total=df.shape[0],
    desc="Extracting Target Classification",
    unit="rows",
):
    target_id = row["chembl_id"]
    if target_id is None or pd.isna(target_id):
        target_classification.append(None)
        continue
    target = new_client.target
    res = target.get(target_id)

    if len(res["target_components"]) != 1:
        print("Warning: More than one component found for target:", target_id)

        # Take only the first component
        res["target_components"] = [res["target_components"][0]]

    comps = res.get("target_components", [])

    if "component_id" in comps[0]:
        component_id = comps[0]["component_id"]

        target_comp = new_client.target_component
        res = target_comp.get(component_id)

        protein_classifications = res.get("protein_classifications", [])

        if "protein_classification_id" in protein_classifications[0]:
            protein_classification_id = protein_classifications[0][
                "protein_classification_id"
            ]

            protein_class = new_client.protein_classification
            class_res = protein_class.get(protein_classification_id)

            if "protein_class_desc" in class_res:
                target_classification.append(class_res["protein_class_desc"])
            else:
                target_classification.append(None)
                print("No protein class description found.")
        else:
            target_classification.append(None)
            print("No protein classification ID found.")
    else:
        target_classification.append(None)
        print("No component ID found.")

df["target_classification"] = target_classification

df.to_csv("results_target_class.csv", index=False)
