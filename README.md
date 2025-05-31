# 🚀 MITOPROJECT

**Author:** Rithwik Nambiar  
**Course:** OMICS Assignment — BIODS  
**ID:** 20244013  

---

## 📌 Objective

1. Process single-cell RNA-seq (scRNA-seq) data for machine learning (ML) analysis.  
2. Compare the expression patterns of nuclear- and mitochondrial-encoded genes in **control vs. LPS-treated** conditions.

**Dataset:** [PMID: 37414801](https://pubmed.ncbi.nlm.nih.gov/37414801), GEO Series: **GSE226488**

---

## 📂 Project Contents

| Folder/File               | Description                                                                 |
|---------------------------|-----------------------------------------------------------------------------|
| `Input/`                  | Contains `lps/` and `control/` subfolders with raw input data              |
| `1_DataPreprocessing.R`   | R script to preprocess data and generate `Mito_expression_data.csv`        |
| `2_ML_Classifiers.py`     | Python script for training ML classifiers using the preprocessed data      |
| `3_GO_Analysis.R`         | Performs Gene Ontology (GO) analysis using top genes from ML output        |
| `Report/`                 | PDF version of the final report and all associated figures                 |
| `Human.MitoCarta3.0/`     | Contains data used in gene selection (OXPHOS and mtDNA-encoded genes)      |

---

## 🧪 How to Run the Project

1. **Preprocess Data**  
   Set your working directory and run `1_DataPreprocessing.R`  
   ➜ Output: `Mito_expression_data.csv`

2. **Train Machine Learning Models**  
   Run `2_ML_Classifiers.py` using the output from Step 1  
   ➜ Output: `top_mito_genes_for_GO.csv`

3. **Perform GO Analysis**  
   Run `3_GO_Analysis.R` using the gene list from Step 2  
   ➜ Output: GO enrichment results

---

## 📋 Notes

- Make sure required R and Python packages are installed before running the scripts.
- For reproducibility, it’s recommended to run each step in a clean environment.
- If you use this project or its structure, please cite appropriately or give credit.

---

