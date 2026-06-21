# Comparation databases (Python)

This Python script compares, within the same study, the total number and percentage of significantly differentially expressed miRNAs detected when contrasting the experimental data against different databases: miRBase, MirGeneDB, and miRCarta.

---

## Input files

- **Directory containing sRNAde results (differential expression matrices)**
  - Selected by the user 
  - The folder for the assignment method must be named: 
	-	de_rcsa 
	-	de_rcadj
  - The matrices must be named: 
	-	`DESeq2_[Condition 1]_vs_[Condition 2].tsv`
	-	`edgeR_[Condition 1]_vs_[Condition 2].tsv`
	-	`limma_[Condition 1]_vs_[Condition 2].tsv`

---

## Requirements

- **Python 3**
- Required libraries:
  - `pandas`
  - `matplotlib`
  - `seaborn`

No additional external dependencies are required.

---

## Script usage

Run the script from the terminal or a Python environment:

```bash
python comparation_data_bases.py -i <input_file> -o <output_file> [-s <study>]
````

Arguments
| Argument | Alternative_Argument | Required | Description | Default |
| --- | --- | --- | --- | --- |
| -i | --input | Yes | Path to the input file | - |
| -o | --output | Yes | Path to the output file | - |
| -s | --study | No | Select the name of the study | "SRP" |

---

## Output

The script generates five output files for every study:

1. **Tab-separated values file (tsv)**
	- Contains the processed data with all the values used in the graphical representations.

2. **Bar plot for significant miRNAs (png)**

3. **Bar plot for percent of significant miRNAs (png)**

4. **Bar plot for significantly differentially expressed miRNAs (png)**

5. **Bar plot for percent of significantly differentially expressed miRNAs (png)**

---

## How the script works

The script is divided into three main stages:

### 1️. File selection and loading

- Searches the user-selected folder for all differential expression matrices and stores them in a list.  

### 2. Generation of the `DataFrame`

- Each file is iterated over individually.
- The number of significant miRNAs and significantly differentially expressed miRNAs, along with their corresponding percentages, are calculated.
- Those values are added to a `DataFrame` containing information for all methodological combinations within a study.

### 3. Graphical representations

- A bar chart summarizing all the information is generated and saved.
	- Significant miRNAs (padj < 0.05)
 	- Significantly differentially expressed miRNAs (padj < 0.05 and abs(log2FoldChange) > 1)
