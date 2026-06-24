# Systematic Computational Identification of Active Drug Metabolites Across the Human Proteome

This repository contains the code used to generate the results reported in the article Larsson et al. (2026) "Computational Identification of Active Drug Metabolites for Human Protein Targets", published open-access in [Molecular Pharmaceutics](https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.6c00258).

Here we conduct a systematic computational analysis to identify drug-metabolite pairs in which the corresponding metabolites retain or increase in bioactivity relative to their parent drugs. 
As part of this work, we used previously curated datasets of drug-metabolite pairs observed in humans ([LAGOM](https://github.com/tsofiac/LAGOM) and [GLORYx](https://github.com/christinadebruynkops/GLORYx))
and related them to bioactivity measurements recorded in [BindingDB](https://www.bindingdb.org) and [ChEMBL](https://www.ebi.ac.uk/chembl/). 
The bioactivity measurements were extracted and rigorusly curated following the procedure described in the paper, with code to reproduce those data extraction and curation steps outlined in this repo.
The dataset was further enhanced with drug-metabolite specific information before an in-depth analysis of the bioactivity trends in the data, reported in the paper.

## Virtual Environment

We recommend to use the [uv](https://github.com/astral-sh/uv) Python package manager, which can be installed following the directions on the linked page. Once installed, you can load the dependencies for the project from ``pyproject.toml``:

```
uv python install 3.10
uv sync
```

## Dataset

Two datasets are used in this project:
* **GLORYx**: 136 drug-metabolite pairs related to the top-selling drugs of 2018 ([gloryx_test_dataset.json](https://github.com/christinadebruynkops/GLORYx/tree/master/datasets/test_dataset)),
* and **LAGOM**:  a rigously curated dataset for drug metabolism, created for [LAGOM](https://github.com/tsofiac/LAGOM)

Obtain the LAGOM dataset, and structure both datasets in the same format by using the following scripts from LAGOM:
* ``load_gloryx.py`` - for GLORYx test dataset,
* ``load_drugbank.py`` - drug-metabolite pairs from Drugbank,
* ``load_metxbiodb.py`` - drug-metabolite pairs from MetXBioDB,
* and ``preprocessing_data`` - to combine and curate the drug-metabolite pairs.

Please see the [LAGOM documentation](https://github.com/tsofiac/LAGOM) for more details.

This results in two files called ``LAGOM_finetune.csv`` and ``gloryx_smiles_clean.csv``.

To obtain the dataset needed for this project, run ``0_get_dataset.py``. This creates a file called ``extended_LAGOM_dataset.csv``.

## Extract Bioactivity Data

For extracting bioactivity data from [ChEMBL](https://www.ebi.ac.uk/chembl/) and [BindingDB](https://www.bindingdb.org/rwd/bind/index.jsp), we will need the dataset file ``extended_LAGOM_dataset.csv``.

* BindingDB: [BindingDB_All.tsv](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp)
* ChEMBL 36: [chembl_36_chemreps.txt](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/)

The following scripts extracts and curates the bioactivity data from the databases BindingDB and ChEMBL:
* ``1_extract_from_bindingdb.py`` - bioactivity data from BindingDB
* ``2_extract_from_chembl.py`` - bioactivity data from ChEMBL

This results in the files ``bindingdb_smiles.csv`` and ``chembl_smiles.csv``.

## Curate Bioactivity Data

For curating the newly extracted biactivity data, we will need the files ``extended_LAGOM_dataset.csv``, ``bindingdb_smiles.csv`` and ``chembl_smiles.csv``.

If you have all those files, please run the notebook: ``3_curate_data.ipynb``. 

It filters, aggregates, and matches the data into a dataset with metabolic reactions and corresponding bioactivity measurements. It also prints some general statistics of the dataset.

This creates two files, named ``matched_reactions_deduplicated.csv`` and ``unique_smiles.csv``.

## Add Drug-Metabolite-Specific Information to Dataset

To add physichemical descriptors and target and reaction information to the dataset, you will need to run three separate scripts:
* ``4_create_SMARTS.py`` - creates reaction templates (SMARTS) from the drug-metabolite pairs
   * Reads: ``matched_reactions_deduplicated.csv``
   * Returns: ``results_SMARTS.py``
* ``4_logd_extraction.py`` - extracts logD7.4 values for each compound predicted by the [ADMETlab 3.0](https://admetlab3.scbdd.com/) model
   * Reads: ``unique_smiles.csv``
   * Returns: ``results_logd.csv``
* ``4_target_class_extraction.py`` - extracts ChEMBL target classifications of each target
   * Reads: ``matched_reactions_deduplicated.csv``
   * Returns: ``results_target_class.csv``

Once collected, all of this new data can be added to the dataset using the script ``5_add_to_data.py``. The following properties are also added to the dataset for each molecule:
* RDKit-calculated logarithmic partition coefficient (cLogP)
* RDKit-calculated topological polar surface area (TPSA)

The filename for the updated dataset with added drug-metabolite-specific information is ``enriched_reactions.csv``.

## Clearance

To predict human hepatocyte CLint for each of the compounds in the drug-metabolite pairs we used an internal AstraZeneca prediction model (not provided).

The following script was run to obtain data in the right format for plotting:
* ``6_clearance.py``
   * Reads: predicted HH CLint values of the ``unique_smiles.csv`` and ``enriched_reactions.csv``
   * Returns: ``results_clearance.csv``

## Analysis

The full analysis of the bioactivity data is summarized in a Python notebook (``7_analysis.ipynb``).

Additional analysis of prodrugs in the dataset was performed using the script called ``8_prodrug_analysis.py``.
 
## Contributors

Sofia Larsson: [@tsofiac](https://github.com/tsofiac)

Rocío Mercado Oropeza: [@rociomer](https://github.com/rociomer)

Susanne Winiwarter

Filip Miljković: [@filipm90](https://github.com/filipm90)
