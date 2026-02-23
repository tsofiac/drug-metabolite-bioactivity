import pandas as pd
import numpy as np
from tqdm import tqdm
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

def extract_unique_smiles(df):

    # Combine SMILES from both drug and metabolite columns
    combined = pd.concat([df["drug"], df["metabolite"]], ignore_index=True)

    # Remove duplicate SMILES strings
    unique_smiles = combined.drop_duplicates()
    
    # Convert to DataFrame
    df = unique_smiles.to_frame(name="smiles")

    return df


def extract_from_bindingdb(df, bindingdb):

    # Convert LAGOM SMILES column to a list
    lagom_smiles = df["smiles"].tolist()

    # Define which columns to extract from the BindingDB file
    columns_to_keep = [
        "Ligand SMILES",
        "Ki (nM)",
        "IC50 (nM)",
        "Kd (nM)",
        "Target Source Organism According to Curator or DataSource",
        "Number of Protein Chains in Target (>1 implies a multichain complex)",
        "UniProt (SwissProt) Primary ID of Target Chain 1",
    ]

    matching_rows = []

    # Set chunk size for reading large BindingDB file in batches
    chunksize = 100_000
    for chunk in pd.read_csv(
        bindingdb, sep="\t", quoting=3, engine="python", chunksize=chunksize
    ):
        # Check if the required SMILES column exists in this chunk
        if "Ligand SMILES" in chunk.columns:
            chunk = chunk[columns_to_keep]
            chunk["Ligand SMILES"] = standardize_smiles_collection(
                chunk["Ligand SMILES"].tolist(), False
            )
            # Filter chunk to keep only rows where SMILES match LAGOM dataset
            filtered_chunk = chunk[chunk["Ligand SMILES"].isin(lagom_smiles)]
            if not filtered_chunk.empty:
                matching_rows.append(filtered_chunk)

    result_df = pd.concat(matching_rows, ignore_index=True)

    print(f"Unique SMILES in LAGOM: {len(set(lagom_smiles))}")
    print(
        f"Unique SMILES in LAGOM with BindingDB data: {result_df['Ligand SMILES'].nunique()}"
    )

    return result_df


def filter_bindingdb_results(df):
    # Filter to keep only single-chain protein targets (exclude complexes)
    df_1 = df[
        (
            df["Number of Protein Chains in Target (>1 implies a multichain complex)"]
            <= 1
        )
    ]
    # Further filter to keep only human targets
    df_2 = df_1[
        df_1["Target Source Organism According to Curator or DataSource"]
        == "Homo sapiens"
    ]

    all_rows = []

    standard_types = ["Ki (nM)", "IC50 (nM)", "Kd (nM)"]

    # Iterate through each row with progress tracking
    for idx, row in tqdm(df_2.iterrows(), total=df_2.shape[0]):
        for standard_type in standard_types:
            value = row.get(standard_type)
            if pd.notna(value):
                if value.startswith(">"):
                    relation = ">"
                    value = value.lstrip(">").strip()
                    value = float(value)
                elif value.startswith("<"):
                    relation = "<"
                    value = value.lstrip("<").strip()
                    value = float(value)
                else:
                    relation = "="
                    try:
                        value = float(value)
                    except ValueError:
                        print(f"ValueError: {value}")
                        continue
                
                # Skip rows without a valid UniProt ID
                if row[
                    "UniProt (SwissProt) Primary ID of Target Chain 1"
                ] is None or pd.isna(
                    row["UniProt (SwissProt) Primary ID of Target Chain 1"]
                ):
                    continue

                clean_type = standard_type.replace(" (nM)", "")
                new_row = {
                    "smiles": row["Ligand SMILES"],
                    "bioactivity": value,
                    "bioactivity_type": clean_type,
                    "bioactivity_relation": relation,
                    "uniprot_id": row[
                        "UniProt (SwissProt) Primary ID of Target Chain 1"
                    ],
                }
                all_rows.append(new_row)

    result_df = pd.DataFrame(all_rows)

    result_df.loc[:, "source"] = "BindingDB"

    print(
        f"Unique SMILES in LAGOM with cleaned BindingDB data: {result_df['smiles'].nunique()}"
    )

    return result_df


if __name__ == "__main__":
    bindingdb = "raw_data/BindingDB_All.tsv"

    lagom_df = pd.read_csv("data/extended_LAGOM_dataset.csv")
    lagom_df = extract_unique_smiles(lagom_df)

    results_df = extract_from_bindingdb(lagom_df, bindingdb)
    results_df.to_csv("data/bindingdb_raw_results.csv", index=False)

    results_df = pd.read_csv("data/bindingdb_raw_results.csv")
    result_df = filter_bindingdb_results(results_df)
    result_df.to_csv("data/bindingdb_smiles.csv", index=False)

    # --- Output extended LAGOM ---
    # Unique SMILES in LAGOM: 5322
    # Unique SMILES in LAGOM with BindingDB data: 1509
    # Faulty entries removed: .10,000
    # Unique SMILES in LAGOM with cleaned BindingDB data: 1169
