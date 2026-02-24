import pandas as pd
from chembl_webresource_client.new_client import new_client


smiles_with_chembl_id = pd.read_csv("data/extended_LAGOM_with_chembl36_ids.csv")
df = pd.read_csv("drug-metabolite_pairs_active.csv")

df_unique = df[["drug", "metabolite"]].drop_duplicates()

length_before = len(df_unique)

# Merge smiles_with_chembl_id to df_unique on drug and metabolite
df_unique = pd.merge(
    df_unique,
    smiles_with_chembl_id[["chembl_ids", "smiles"]],
    left_on="drug",
    right_on="smiles",
    how="left",
)
df_unique = df_unique.rename(columns={"chembl_ids": "drug_chembl_ids"})
df_unique = df_unique.drop(columns=["smiles"])
df_unique = pd.merge(
    df_unique,
    smiles_with_chembl_id[["chembl_ids", "smiles"]],
    left_on="metabolite",
    right_on="smiles",
    how="left",
)
df_unique = df_unique.rename(columns={"chembl_ids": "metabolite_chembl_ids"})
df_unique = df_unique.drop(columns=["smiles"])

df_unique.to_csv("prodrugs.csv", index=False)


# Create a molecule client
molecule = new_client.molecule


for index, row in df_unique.iterrows():
    drug_smiles = row["drug"]
    metabolite_smiles = row["metabolite"]
    drug_chembl_ids = row["drug_chembl_ids"]

    # Interpret chembl_ids as lists
    drug_chembl_ids = eval(drug_chembl_ids)

    # Initialize as not a prodrug
    df_unique.at[index, "is_prodrug"] = False

    id_list = []
    for drug_chembl_id in drug_chembl_ids:

        # Get the molecule data
        mol_data = molecule.get(drug_chembl_id)

        if mol_data and "prodrug" in mol_data and mol_data["prodrug"] == 1 and "molecule_hierarchy" in mol_data:
            hierarchy = mol_data["molecule_hierarchy"]
            active_chembl_id = hierarchy.get("active_chembl_id")

            # Check if metabolite_chembl_ids contains parent_chembl_id or active_chembl_id
            if active_chembl_id in eval(row["metabolite_chembl_ids"]):
                df_unique.at[index, "is_prodrug"] = True
                id_list.append(drug_chembl_id)
                print(f"For prodrug ChEMBL ID: {drug_chembl_id}")
                print(f"    Active ChEMBL ID: {active_chembl_id}")
                print(
                    f"    Metabolite ChEMBL IDs: {eval(row['metabolite_chembl_ids'])}"
                )

    df_unique.at[index, "prodrug_chembl_ids"] = str(id_list)

df_prodrugs = df_unique[df_unique["is_prodrug"]]
df_prodrugs.to_csv("prodrugs.csv", index=False)


mechanism = new_client.mechanism
for index, row in df_prodrugs.iterrows():
    prodrug_chembl_ids = eval(row["prodrug_chembl_ids"])
    metabolite_chembl_ids = eval(row["metabolite_chembl_ids"])

    moa_drug = []
    moa_metabolite = []
    moa_chembl_ids_drug = []
    moa_chembl_ids_metabolite = []
    for prodrug_chembl_id in prodrug_chembl_ids:
        mechanisms_prodrug = mechanism.filter(molecule_chembl_id=prodrug_chembl_id)
        mechanisms_metabolite = mechanism.filter(
            molecule_chembl_id__in=active_chembl_id
        )

        for mech_d in mechanisms_prodrug:
            if mech_d:
                # Get mechanism of action description
                action_d = mech_d.get("mechanism_of_action")
                print(f"        Mechanism of action for drug: {action_d}")
                if action_d not in moa_drug:
                    moa_drug.append(action_d)
                # Get target_chembl_id
                target_chembl_id = mech_d.get("target_chembl_id")
                if target_chembl_id not in moa_chembl_ids_drug:
                    moa_chembl_ids_drug.append(target_chembl_id)
        for mech_m in mechanisms_metabolite:
            if mech_m:
                # Get mechanism of action description
                action_m = mech_m.get("mechanism_of_action")
                print(f"        Mechanism of action for metabolite: {action_m}")
                if action_m not in moa_metabolite:
                    moa_metabolite.append(action_m)
                # Get target_chembl_id
                target_chembl_id = mech_m.get("target_chembl_id")
                if target_chembl_id not in moa_chembl_ids_metabolite:
                    moa_chembl_ids_metabolite.append(target_chembl_id)

    df_prodrugs.at[index, "moa_drug"] = str(moa_drug)
    df_prodrugs.at[index, "moa_metabolite"] = str(moa_metabolite)
    df_prodrugs.at[index, "moa_chembl_ids_drug"] = str(moa_chembl_ids_drug)
    df_prodrugs.at[index, "moa_chembl_ids_metabolite"] = str(moa_chembl_ids_metabolite)

df_prodrugs.to_csv("prodrugs.csv", index=False)

df_active = pd.read_csv("enriched_reactions_active.csv")

for index, row in df_prodrugs.iterrows():
    drug_smiles = row["drug"]
    metabolite_smiles = row["metabolite"]

    # Find rows in df_active that match drug_smiles and metabolite_smiles
    matching_rows = df_active[
        (df_active["drug"] == drug_smiles)
        & (df_active["metabolite"] == metabolite_smiles)
    ]
    target_ids = []
    # Append the unique target_ids from matching_rows to target_ids list
    for _, match_row in matching_rows.iterrows():
        target_id = match_row["chembl_id"]
        if target_id not in target_ids:
            target_ids.append(target_id)

    print(target_ids)

    df_prodrugs.at[index, "target_ids"] = target_ids

df_prodrugs.to_csv("prodrugs.csv", index=False)

df_prodrugs = pd.read_csv("prodrugs.csv")

target_relation = new_client.target_relation
count_true = 0
count_all = 0

for index, row in df_prodrugs.iterrows():
    moa_chembl_ids_drug = eval(row["moa_chembl_ids_drug"])
    moa_chembl_ids_metabolite = eval(row["moa_chembl_ids_metabolite"])
    target_ids = eval(row["target_ids"])

    # Check for direct overlap first, then check for target relations if no direct overlap is found
    overlap_results = []
    for target_id in target_ids:
        overlap = False
        count_all += 1
        if target_id in moa_chembl_ids_drug:
            overlap = True
            count_true += 1
        if target_id in moa_chembl_ids_metabolite:
            overlap = True
            count_true += 1

        if overlap:
            print(
                f"    Overlap found for target ID {target_id} with MOA ChEMBL IDs: {moa_chembl_ids_drug + moa_chembl_ids_metabolite}"
            )
        else:
            target_relations = target_relation.filter(
                target_chembl_id=moa_chembl_ids_drug + moa_chembl_ids_metabolite,
                relation_type="SAME_TARGET_CLASS",
                related_target_chembl_id=target_id,
            )
            if target_relations:
                print(
                    f"    Target relation found for target ID {target_id} with MOA ChEMBL IDs: {moa_chembl_ids_drug + moa_chembl_ids_metabolite}"
                )
                overlap = True
                count_true += 1
            else:
                print(
                    f"    No overlap or target relation found for target ID {target_id} with MOA ChEMBL IDs: {moa_chembl_ids_drug + moa_chembl_ids_metabolite}"
                )

        overlap_results.append(overlap)

    df_prodrugs.at[index, "moa_target_overlap"] = str(overlap_results)

df_prodrugs.to_csv("prodrugs.csv", index=False)


print("Summary:")
print(
    f"Total number of unique prodrugs identified: {len(df_prodrugs)} out of {length_before} unique drug-metabolite pairs."
)
print(f"Total overlaps found: {count_true} out of {count_all} reaction-target pairs")
