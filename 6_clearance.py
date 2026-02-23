import pandas as pd

df_predicted = pd.read_csv(r".../PIP_Results_2026-01-26_10-06-21.csv")

df = pd.read_csv("enriched_reactions.csv")
df = df[["drug", "metabolite"]]
df = df.drop_duplicates()

df_predicted = df_predicted[["smiles", "HH CLint_HH logCLint"]]
df_predicted = df_predicted.rename(columns={"HH CLint_HH logCLint": "hh"})
df = df.merge(
    df_predicted.rename(columns={"smiles": "drug"}),
    on="drug",
    how="left",
)
df = df.rename(columns={"hh": "hh_drug"})
df = df.merge(
    df_predicted.rename(columns={"smiles": "metabolite"}), on=["metabolite"], how="left"
)
df = df.rename(columns={"hh": "hh_metabolite"})
df = df[["drug", "metabolite", "hh_drug", "hh_metabolite"]]

df.to_csv(r".../results_clearance.csv", index=False)
