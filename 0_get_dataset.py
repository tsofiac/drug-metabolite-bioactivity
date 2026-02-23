import pandas as pd


def combine_datasets(LAGOM, GLORYx):

    # Combine all dataframes
    combined_df = pd.concat([LAGOM, GLORYx], ignore_index=True)

    # Find duplicate pairs and show which sources they came from
    duplicates = combined_df[combined_df.duplicated(subset=['drug', 'metabolite'], keep=False)]
    if not duplicates.empty:
        print("Found duplicate drug-metabolite pairs:")
        # Group by drug-metabolite pairs and show their sources
        for (drug, metabolite), group in duplicates.groupby(['drug', 'metabolite']):
            print(f"\nPair: {drug},{metabolite}")
            print("Found in:", group['file'].tolist())
    else:
        print("No duplicate drug-metabolite pairs found")

    # Keep only unique pairs in the combined dataframe
    combined_df = combined_df.drop_duplicates(subset=['drug', 'metabolite'])
    
    combined_df = combined_df[['drug', 'metabolite', 'source']].reset_index(drop=True)

    return combined_df

    

if __name__ == "__main__":

    # CREATE LAGOM DATASET
    # ADD GLORYx DATASET 

    # Before running this code, make sure you have downloaded the LAGOM_finetune.csv file as well as gloryx_smiles_clean.csv
    # from the LAGOM repository: github.com/tsofiac/LAGOM
    # Use the load_gloryx.py script to create the gloryx_smiles_clean.csv file from the gloryx_test_dataset.json file
    # Use load_drugbank.py and load_metxbiodb.py and preprocessing_data.py and standardize_smiles.py to create the LAGOM_finetune.csv file from DrugBank and MetXBioDB datasets


    LAGOM = pd.read_csv("LAGOM/LAGOM_finetune.csv", sep="\t")
    LAGOM = LAGOM.rename(columns={"reactants": "drug", "products": "metabolite"})
    print(f"LAGOM dataset size: {len(LAGOM)}")
    GLORYx = pd.read_csv("LAGOM/gloryx_smiles_clean.csv")
    GLORYx = GLORYx.rename(columns={"parent_smiles": "drug", "child_smiles": "metabolite"})
    print(f"GLORYx dataset size: {len(GLORYx)}")

    df = combine_datasets(LAGOM, GLORYx)
    print(f"Extended LAGOM dataset size: {len(df)}")
    df.to_csv("data/extended_LAGOM_dataset.csv", index=False)

