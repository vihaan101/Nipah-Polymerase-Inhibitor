# Nipah Virus Inhibitor Discovery Kit
**Clean, Simplified, and Reproducible.**

This repository contains the verified pipeline for discovering inhibitors for the Nipah Virus Polymerase (L-Protein), specifically targeting the W730 allosteric pocket to ensure resilience against the common W730A escape mutation.

## Directory Structure
- `Experiment_Notebook.ipynb`: **Start Here.** The interactive lab notebook that explains every step in plain English.
- `pipeline.py`: The tool for you to test your own drug candidates.
- `verify.py`: The "Certificate of Authenticity" script to verify our findings (BMS-986205 vs ERDRP-0519).
- `data/`: Contains the precise 9KNZ receptor structures (Wild-Type and W730A Mutant) and the ligand files.
- `bin/`: Contains the AutoDock Vina executable (macOS Silicon/Intel compatible).
- `HOW_I_BUILT_THIS.md`: The "Director's Cut" of the project—chemistry choices, failures, and fixes.
- `Final_Report.md`: The formal scientific paper explaining the biology, methodology, and results.

## 1. Quick Start
**Recommended:** Open `Experiment_Notebook.ipynb` in Jupyter or VS Code. It runs the entire experiment interactively and plots the graphs for you.

If you prefer the command line:
If you just want to run the pipeline immediately, the `data/` folder already contains the pre-computed receptor files.

```bash
python3 verify.py
```

## 2. The "From Scratch" Route (Pure Reproducibility)
If you want to verify the entire data preparation pipeline starting from the raw PDB files (e.g., for an audit):

1.  **Install Prerequisites:**
    -   **OpenBabel**: You must have `obabel` in your path.
        -   Mac: `brew install open-babel`
        -   Linux: `sudo apt-get install openbabel`
    -   **Python Deps**: `pip install biopython gemmi numpy pandas meeko rdkit`

2.  **Run the Setup:**
    ```bash
    python3 setup.py
    ```
    This will:
    -   Download the raw `9KNZ.cif` from the Protein Data Bank.
    -   Clean the protein structure computationally.
    -   Generate the W730A mutant structure *in silico*.
    -   Convert everything to the PDBQT physics format used by Vina.

## 3. How to Discover Your Own Inhibitor
If you have a molecule (in PDBQT format) that you want to test against the Nipah virus:

```bash
python3 pipeline.py path/to/your_ligand.pdbqt
```

**What this does:**
1.  Docks your molecule into the Wild-Type Nipah Polymerase.
2.  Docks it into the Resistant Mutant (W730A).
3.  Calculates the "Resistance Delta" (how much potency is lost).
4.  Gives you a pass/fail verdict based on **Potency** (<-6.5 kcal/mol) and **Resilience** (<0.5 kcal/mol loss).

## 4. How to Verify Our Results
To reproduce the scientific claims that **BMS-986205** is resistant while **ERDRP-0519** is susceptible:

```bash
python3 verify.py
```

**What this does:**
-   Runs a rigorous statistical simulation using 5 unique random seeds.
-   Performs a head-to-head comparison of BMS-986205 vs ERDRP-0519.
-   Outputs a final "Certificate of Reproducibility" proving the resistance profile.

## 5. Understanding the Science (The "Why")
This code is just the final product. But the real science happened in the failures. To read about the "Vacuum Hole Fallacy", the "Crystal Water Trap", and other rabbit holes I went down (and how I climbed out), check out:

👉 **[HOW_I_BUILT_THIS.md](HOW_I_BUILT_THIS.md)**

It's the personal story of how this project was actually built.

## 6. The Formal Paper
If you want the standard scientific format (Abstract, Methods, Results, References), see:

👉 **[Final_Report.md](Final_Report.md)**

## Requirements
-   **Python 3.9+**
-   **macOS** (This kit includes a macOS binary for Vina. For Linux, replace `bin/vina` with a Linux binary).
-   Standard libraries: `sys`, `subprocess`, `statistics`, `pathlib`.
