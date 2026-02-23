import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen


def add_lagom_source(df, lagom_df):

    lagom_df = lagom_df[["drug", "metabolite", "source"]]
    lagom_df = lagom_df.rename(columns={"source": "source_lagom"})
    df = df.merge(
        lagom_df,
        left_on=["drug", "metabolite"],
        right_on=["drug", "metabolite"],
        how="left",
    )
    print(df["source_lagom"].value_counts().to_string())

    return df


def add_molecular_properties(df):
    smiles = df["drug"].tolist() + df["metabolite"].tolist()
    mols = [Chem.MolFromSmiles(x) for x in smiles]

    logp = [round(Crippen.MolLogP(x), 2) for x in mols]
    tpsa = [round(Descriptors.TPSA(x), 2) for x in mols]

    df = df.copy()
    df.loc[:, "logp_drug"] = logp[: len(df)]
    df.loc[:, "logp_metabolite"] = logp[len(df) :]
    df.loc[:, "tpsa_drug"] = tpsa[: len(df)]
    df.loc[:, "tpsa_metabolite"] = tpsa[len(df) :]

    return df


def add_logd_to_reactions(logd_results, df):

    df1 = df
    df2 = pd.read_csv(logd_results)

    df1 = df1.merge(df2, left_on="drug", right_on="smiles", how="left")
    df1 = df1.rename(columns={"logD7.4": "logd_drug"})
    df1 = df1.drop(columns=["smiles"])

    df1 = df1.merge(df2, left_on="metabolite", right_on="smiles", how="left")
    df1 = df1.rename(columns={"logD7.4": "logd_metabolite"})
    df1 = df1.drop(columns=["smiles"])

    return df1


def add_target_class_to_reactions(target_class_results, df):

    df1 = df
    df2 = pd.read_csv(target_class_results)

    final_df = pd.merge(df1, df2, on="uniprot_id", how="left")
    final_df["target_classification"] = final_df["target_classification"].fillna(
        "missing"
    )

    return final_df


def add_reaction_template(templates_file, df):

    df_templates = pd.read_csv(templates_file)
    df_templates = df_templates[["template r=0", "drug", "metabolite"]]
    print(f"Number of unique templates: {df_templates['template r=0'].nunique()}")
    print(
        f"Number of unique drug-metabolite pairs: {df_templates[['drug', 'metabolite']].shape[0]}"
    )
    df = df.merge(df_templates, on=["drug", "metabolite"], how="left")

    return df


if __name__ == "__main__":
    df = pd.read_csv("inbetween_data/matched_reactions_deduplicated.csv")
    print(f"Initial number of reactions: {df.shape[0]}")

    lagom_df = pd.read_csv("data/extended_LAGOM_dataset.csv")
    df = add_lagom_source(df, lagom_df)

    df = add_molecular_properties(df)

    logd_results = "results_logd.csv"
    df = add_logd_to_reactions(logd_results, df)

    target_class_results = "results_target_class.csv"
    df = add_target_class_to_reactions(target_class_results, df)

    templates_file = "results_smarts.csv"
    df = add_reaction_template(templates_file, df)

    print(f"Final number of reactions: {df.shape[0]}")
    df.to_csv("enriched_reactions.csv", index=False)
