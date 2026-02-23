# import os
import ast
from time import sleep
import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm
import numpy as np
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from chembl_structure_pipeline.standardizer import standardize_mol


def standardize_molecule(smiles, isomericSmiles=True):
    """
    Standardization of SMILES by standardiser;
    """

    try:
        rm = SaltRemover()
        mol = Chem.MolFromSmiles(smiles)
        mol = rm.StripMol(mol, dontRemoveEverything=True)
        if len(Chem.MolToSmiles(mol).split(".")) > 1:
            salt = Chem.MolToSmiles(mol)
            frag_dict = {len(k): k for k in salt.split(".")}
            max_frag = frag_dict[max(list(frag_dict.keys()))]
            mol = Chem.MolFromSmiles(max_frag)
        mol = standardize_mol(mol)
        clean_smiles = Chem.MolToSmiles(
            mol, isomericSmiles=isomericSmiles
        )  # removing (False) or keeping (True) stereochemistry
        if Chem.MolFromSmiles(clean_smiles) is None:
            clean_smiles = None
    except Exception:
        clean_smiles = np.nan
    return clean_smiles


def standardize_smiles_collection(smiles_list, isomericSmiles=True):
    """
    Standardization of SMILES collection by standardiser;
    """

    lookup = {}
    return_smiles = []
    for smiles in smiles_list:
        if smiles in lookup:
            standardized_smiles = lookup[smiles]
        else:
            standardized_smiles = standardize_molecule(
                smiles, isomericSmiles=isomericSmiles
            )
            lookup[smiles] = standardized_smiles
        return_smiles.append(standardized_smiles)

    return return_smiles


def standardize_chembl_smiles(input_file, output_file):

    # Read ChEMBL data
    chembl_df = pd.read_csv(input_file, sep="\t", usecols=["chembl_id", "canonical_smiles"])

    smiles = chembl_df["canonical_smiles"].tolist()
    standardized = []

    # Process SMILES in batches
    batch_size = 10000
    for i in tqdm(range(0, len(smiles), batch_size), desc="Standardizing ChEMBL SMILES", unit="molecules"):
        chunk = smiles[i:i+batch_size]
        standardized.extend(standardize_smiles_collection(chunk, isomericSmiles=False))

    # Add standardized SMILES as a new column
    chembl_df["standardized_smiles"] = standardized
    chembl_df.to_csv(output_file, sep="\t", index=False)    


def extract_unique_smiles(df):

    # Combine SMILES from both drug and metabolite columns
    combined = pd.concat([df["drug"], df["metabolite"]], ignore_index=True)

    # Remove duplicate SMILES strings
    unique_smiles = combined.drop_duplicates()
    
    # Convert to DataFrame
    df = unique_smiles.to_frame(name="smiles")

    return df

def match_chembl_ids_to_lagom(chembl_std, lagom_df):

    # Load standardized ChEMBL data
    chembl_df = pd.read_csv(chembl_std, sep="\t")

    # Match LAGOM SMILES with ChEMBL IDs
    merged_df = pd.merge(
        lagom_df,
        chembl_df[["chembl_id", "standardized_smiles"]],
        left_on="smiles",
        right_on="standardized_smiles",
        how="left"
    )

    # Group by SMILES and aggregate ChEMBL IDs into lists
    grouped = (
        merged_df.groupby("smiles")["chembl_id"]
        .apply(lambda x: list(x.dropna().unique()))
        .reset_index()
    )

    # Merge grouped ChEMBL IDs back to original LAGOM
    final_df = lagom_df.merge(grouped, on="smiles", how="left")
    final_df = final_df.rename(columns={"chembl_id": "chembl_ids"})
    
    # Keep only those with at least one ChEMBL ID
    final_df = final_df[final_df["chembl_ids"].map(len) > 0]
    print(f"Unique SMILES in LAGOM with ChEMBL IDs: {len(final_df)}")

    return final_df

def get_filtered_bioactivity(molecule_id):

    # Initialize ChEMBL API clients
    activity = new_client.activity
    assay = new_client.assay
    target = new_client.target

    # Query bioactivity data with initial filters
    activities = activity.filter(
        molecule_chembl_id=molecule_id, 
        standard_type__in=["IC50", "Ki", "Kd"], 
        standard_units__in=["nM", "pM", "uM", "mM", "M"]).only([
        "molecule_chembl_id",
        "assay_chembl_id",
        "target_chembl_id",
        "standard_type",
        "standard_value",
        "standard_units",
        "standard_relation",
        "activity_comment"
    ])

    filtered = []

    # Apply additional filters for each activity record
    for act in activities:
        target_id = act.get("target_chembl_id")
        assay_id = act.get("assay_chembl_id")

        # Filter 1: Target must be a single human protein
        t_res = list(target.filter(
            target_chembl_id=target_id,
            target_type="SINGLE PROTEIN",
            organism="Homo sapiens"
        ))
        if not t_res:
            continue
            
        # Filter 2: Assay must be direct (D), high confidence (9), and cell-free
        a_res = list(assay.filter(
            assay_chembl_id=assay_id,
            target_chembl_id=target_id,
            relationship_type="D",
            confidence_score=9,
            assay_cell_type=None
        ))
        if not a_res:
            continue

        filtered.append(act)

    return filtered

def extract_bioactivity_data(df):

    # Configuration for fresh start
    start_index = 0
    result_df = pd.DataFrame()

    # Uncomment these lines to resume from a previous run if API fails
    # start_index = 2307
    # result_df = pd.read_csv("data/chembl_smiles_tmp.csv")

    # Iterate through each molecule in the dataset
    for idx, row in tqdm(df.iloc[start_index:].iterrows(), total=len(df), desc="Extracting bioactivity data", unit="molecules"):
        chembl_ids = row['chembl_ids']
        chembl_ids = ast.literal_eval(chembl_ids)

        rows = []
        # Process each ChEMBL ID for the current molecule
        for chembl_id in chembl_ids:

            # Retry logic to handle API failures
            max_retries = 5
            for attempt in range(max_retries):
                try:
                    bioactivities = get_filtered_bioactivity(chembl_id)
                    print(f"Found {len(bioactivities)} bioactivities for molecule {chembl_id}")
                    break
                except Exception as e:
                    wait_time = 5 * (attempt + 1)  # 5, 10, 15 seconds
                    if attempt < max_retries - 1:
                        print(f"Error (attempt {attempt + 1}/{max_retries}). Retrying in {wait_time} seconds...")
                        sleep(wait_time)
                    else:
                        print(f"Failed after {max_retries} attempts for molecule {chembl_id}: {e}")
                        exit(1)

            # Process each bioactivity record
            for record in bioactivities:

                # Skip records with missing values or unrecognized units, and convert all values to nM
                if (record['standard_units'] is None) | (record['standard_value'] is None):
                    continue
                elif record['standard_units'] == 'M':
                    print("Molar unit found!")
                    value = float(record['standard_value']) * 1e9  # Convert M to nM
                elif record['standard_units'] == 'mM':     
                    print("Milli molar unit found!")
                    value = float(record['standard_value']) * 1e6  # Convert mM to nM
                elif record['standard_units'] == 'uM':
                    print("Micro molar unit found!")
                    value = float(record['standard_value']) * 1e3  # Convert uM to nM
                elif record['standard_units'] == 'nM':
                    value = float(record['standard_value'])  # nM, no conversion needed
                elif record['standard_units'] == 'pM':
                    print("Pico molar unit found!")
                    value = float(record['standard_value']) * 1e-3  # Convert pM to nM
                else:
                    continue  # Skip if units are unrecognized

                # Create a new row with standardized bioactivity data
                new_row = {
                    'smiles': row['smiles'],
                    'chembl_id': chembl_id,
                    'bioactivity': value,
                    'bioactivity_type': record['standard_type'],
                    'bioactivity_relation': record['standard_relation'],
                    'target_id': record['target_chembl_id'],
                    'activity_comment': record['activity_comment'],
                }

                rows.append(new_row)

        result_df = pd.concat([result_df, pd.DataFrame(rows)], ignore_index=True)
        result_df.to_csv("data/chembl_smiles_tmp.csv", index=False)
        print("Saved temporary results from index:", idx)

    # Map ChEMBL target IDs to UniProt accession numbers
    uniprot_ids = []
    for _, row in tqdm(result_df.iterrows(), total=result_df.shape[0], desc="Extracting UniProt IDs", unit="rows"):
        target_id = row['target_id']
        target = new_client.target
        res = target.get(target_id)

        if len(res["target_components"]) != 1:
            print("Warning: More than one component found for target:", target_id)

        comps = res.get("target_components", [])

        # Extract UniProt accession if available
        if comps and "accession" in comps[0]:
            uniprot_ids.append(comps[0]["accession"])
        else:
            uniprot_ids.append(None)

    result_df['uniprot_id'] = uniprot_ids

    result_df.loc[:, 'source'] = 'ChEMBL'

    print("Unique SMILES in LAGOM with bioactivity data:", result_df['smiles'].nunique())

    return result_df

    

    
if __name__ == "__main__":

    standardize_chembl = False
    if standardize_chembl:
        input_file = "raw_data/chembl_36_chemreps.txt"
        output_file = "raw_data/chembl_36_chemreps_standardized.txt"
        standardize_chembl_smiles(input_file, output_file)

    chembl_std = "raw_data/chembl_36_chemreps_standardized.txt"
    lagom_df = pd.read_csv("data/extended_LAGOM_dataset.csv")
    lagom_df = extract_unique_smiles(lagom_df)
    print(f"Unique SMILES in LAGOM: {len(lagom_df)}")
    lagom_df = match_chembl_ids_to_lagom(chembl_std, lagom_df)
    lagom_df.to_csv("data/extended_LAGOM_with_chembl36_ids.csv", index=False)

    lagom_df = pd.read_csv("data/extended_LAGOM_with_chembl36_ids.csv")
    result_df = extract_bioactivity_data(lagom_df)
    result_df.to_csv("data/chembl_smiles_36.csv", index=False)


    # --- Output extended LAGOM ---
    # Unique SMILES in LAGOM: 5322
    # Unique SMILES in LAGOM with ChEMBL IDs: 2361
    # Unique SMILES in LAGOM with bioactivity data: 1035
