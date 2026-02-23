import pandas as pd
from tqdm import tqdm
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

from rxnmapper import RXNMapper
from rxnutils.chem.reaction import ChemicalReaction


def move_atom_maps_to_notes(m):
    for at in m.GetAtoms():
        if at.GetAtomMapNum():
            at.SetProp("atomNote", str(at.GetAtomMapNum()))


def draw_chemical_reaction(smiles, highlightByReactant=False, font_scale=1.5):
    rxn = rdChemReactions.ReactionFromSmarts(smiles, useSmiles=True)
    trxn = rdChemReactions.ChemicalReaction(rxn)
    # Move atom maps to be annotations:
    for m in trxn.GetReactants():
        move_atom_maps_to_notes(m)
    for m in trxn.GetProducts():
        move_atom_maps_to_notes(m)
    d2d = rdMolDraw2D.MolDraw2DSVG(800, 300)
    d2d.drawOptions().annotationFontScale = font_scale
    d2d.DrawReaction(trxn, highlightByReactant=highlightByReactant)

    d2d.FinishDrawing()

    return d2d.GetDrawingText()


def plot_and_print(mapped_rxn, conf, template, i):
    # Only generate SVG files for the first 10 reactions
    if i < 10:
        svg_data = draw_chemical_reaction(mapped_rxn, highlightByReactant=True)
        filename = f"reaction_template_{i + 1}.svg"
        with open(filename, "w") as f:
            f.write(svg_data)

    print(f"\n=== Reaction {i + 1} ===")
    print(f"Mapped reaction: {mapped_rxn}")
    print(f"Confidence: {conf:.2f}")
    print(f"Template: {template['reaction_smarts']}")


def atom_mapping(df):
    # Create reaction SMILES by combining drug and metabolite with ">>" separator
    df["reaction"] = df["drug"] + ">>" + df["metabolite"]
    rxn_mapper = RXNMapper()
    rxns = df["reaction"].tolist()
    # Generate atom-to-atom mappings using attention-guided algorithm
    results = rxn_mapper.get_attention_guided_atom_maps(rxns)
    mapped_rxns = [result["mapped_rxn"] for result in results]
    confidences = [result["confidence"] for result in results]

    df["mapped_reaction"] = mapped_rxns
    df["confidence"] = confidences

    return df


def extract_template(df):
    # Convert mapped reactions to list for processing
    mapped_rxns = df["mapped_reaction"].tolist()
    confidences = df["confidence"].tolist()

    templates_0 = []
    # Iterate through mapped reactions
    for i, (mapped_rxn, conf) in tqdm(
        enumerate(zip(mapped_rxns, confidences)), total=len(mapped_rxns)
    ):
        rxn = ChemicalReaction(mapped_rxn)
        # Generate reaction template with radius 0 (only reacting atoms, no neighbors)
        rxn.generate_reaction_template(radius=0)
        template_0 = rxn.canonical_template.smarts
        templates_0.append(template_0)

        # Option to plot and print first ten reactions
        # plot_and_print(mapped_rxn, conf, template, i)

    reactions = pd.DataFrame(
        {
            "drug": df["drug"],
            "metabolite": df["metabolite"],
            "reaction": df["reaction"],
            "mapped_reaction": df["mapped_reaction"],
            "confidence": df["confidence"],
            "template r=0": templates_0,
        }
    )

    return reactions


if __name__ == "__main__":
    df = pd.read_csv("inbetween_data/matched_reactions_deduplicated.csv")
    df = df[["drug", "metabolite"]].drop_duplicates().reset_index(drop=True)
    df = atom_mapping(df)  # approx 2 min
    reactions = extract_template(df)  # approx 2h
    reactions.to_csv("results_smarts.csv", index=False)
