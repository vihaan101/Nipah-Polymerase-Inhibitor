
import sys
import subprocess
import json
import re
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results/verification"
SCRIPTS_DIR = BASE_DIR / "scripts"
VINA_PATH = SCRIPTS_DIR / "vina"

# Files
WT_RIGID = DATA_DIR / "9KNZ_clean_v3_hybrid.pdbqt"
MUT_RIGID = DATA_DIR / "9KNZ_W730A_v2.pdbqt_rigid.pdbqt"
LIGAND_MULTI = RESULTS_DIR / "BMS-986205_flex_wt.pdbqt"
LIGAND_SINGLE = RESULTS_DIR / "BMS-986205_pure.pdbqt"
OUTPUT_WT_FILE = RESULTS_DIR / "BMS-986205_rigid_wt.pdbqt"
OUTPUT_MUT_FILE = RESULTS_DIR / "BMS-986205_rigid_mut.pdbqt"

# Config
with open(BASE_DIR / "config/docking_box.json") as f:
    BOX_CONFIG = json.load(f)

def extract_pure_ligand(multimodel_path, output_path):
    print(f"Extracting PURE ligand from {multimodel_path}...")
    with open(multimodel_path, 'r') as f_in:
        content = f_in.readlines()
        
    with open(output_path, 'w') as f_out:
        capture = False
        for line in content:
            if line.startswith("MODEL 1"):
                capture = True
                continue # Skip MODEL tag
            
            if capture:
                # STOP conditions
                if line.startswith("ENDMDL"): break
                if line.startswith("BEGIN_RES"): break 
                
                f_out.write(line)

def parse_output(pdbqt_file):
    best_aff = 999.0
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if "REMARK VINA RESULT" in line:
                val = float(line.split()[3])
                if val < best_aff: best_aff = val
    return best_aff if best_aff != 999.0 else None

def run_docking(receptor, ligand, output):
    print(f"Docking into {receptor.name}...")
    cmd = [
        str(VINA_PATH),
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", str(BOX_CONFIG["center_x"]),
        "--center_y", str(BOX_CONFIG["center_y"]),
        "--center_z", str(BOX_CONFIG["center_z"]),
        "--size_x", str(BOX_CONFIG["size_x"]),
        "--size_y", str(BOX_CONFIG["size_y"]),
        "--size_z", str(BOX_CONFIG["size_z"]),
        "--exhaustiveness", "16", # High precision for final check
        "--scoring", "vinardo",
        "--seed", "42",
        "--out", str(output)
    ]
    print("Running command: " + " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Vina Output:\n{e.stdout.decode()}")
        print(f"Vina Error:\n{e.stderr.decode()}")
        sys.exit(1)

def main():
    print("=== Determining Resistance Profile for BMS-986205 (Rigid vs Rigid Validation) ===")
    
    # 1. Extract pure ligand
    extract_pure_ligand(LIGAND_MULTI, LIGAND_SINGLE)
    
    # 2. Dock WT (Rigid)
    run_docking(WT_RIGID, LIGAND_SINGLE, OUTPUT_WT_FILE)
    
    # 3. Dock Mutant (Rigid)
    run_docking(MUT_RIGID, LIGAND_SINGLE, OUTPUT_MUT_FILE)

    # 4. Compare
    wt_score = parse_output(OUTPUT_WT_FILE)
    mut_score = parse_output(OUTPUT_MUT_FILE)
    
    print(f"WT Score (Rigid):     {wt_score} kcal/mol")
    print(f"Mutant Score (Rigid): {mut_score} kcal/mol")
    
    delta = mut_score - wt_score
    print(f"Delta Affinity:       {delta:.2f} kcal/mol")
    
    if delta < 0.5:
        print("VERDICT: RESISTANT (Passes Energy Check)")
        print("Confirms that the 7A distance translates to energy independence.")
    else:
        print("VERDICT: SUSCEPTIBLE (Fails Energy Check)")

if __name__ == "__main__":
    main()
