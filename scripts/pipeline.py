
import sys
import subprocess
import argparse
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import distance
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

# =============================================================================
# CONFIGURATION
# =============================================================================
# The grid box covering the allosteric pocket near W730
BOX_CONFIG = {
    "center_x": 133.301,
    "center_y": 137.79,
    "center_z": 150.667,
    "size_x": 22,
    "size_y": 22,
    "size_z": 22,
    "exhaustiveness": 16, # High precision
    "num_modes": 1
}

# Paths (Relative to this script)
# Paths (Relative to this script)
# Script is in root/scripts/, so BASE_DIR (root) is parent.parent
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
BIN_DIR = BASE_DIR / "bin"
VINA_EXEC = BIN_DIR / "vina"
WT_RECEPTOR = DATA_DIR / "receptor_wt.pdbqt"
MUT_RECEPTOR = DATA_DIR / "receptor_mut.pdbqt"
WT_PDB = DATA_DIR / "receptor_wt.pdb" # Required for ghost check

# =============================================================================
# FUNCTIONS
# =============================================================================

def check_ghost_clash(ligand_pdbqt, wt_pdb_path=WT_PDB):
    """
    Checks if the ligand occupies the space where the W730 sidechain used to be.
    This prevents the 'Vacuum Hole Fallacy'.
    """
    if not wt_pdb_path.exists():
        return "SKIP" # Cannot check without WT PDB

    # 1. Get Ghost Atoms (W730 Sidechain)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("WT", str(wt_pdb_path))
    ghost_coords = []
    try:
        # Assuming Chain A, Residue 730 (Standard numbering for Nipah)
        # We search chain A or first available chain
        chain = structure[0]['A'] if 'A' in structure[0] else list(structure[0])[0]
        # Find 730
        res = None
        for r in chain:
            if r.id[1] == 730: 
                res = r
                break
        
        if res and res.resname == "TRP":
            for atom in res:
                if atom.name not in ["N", "CA", "C", "O"]: # Sidechain only
                    ghost_coords.append(atom.get_coord())
    except Exception:
        pass

    if not ghost_coords:
        return "SKIP"

    # 2. Get Ligand Atoms from PDBQT
    lig_coords = []
    with open(ligand_pdbqt, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    lig_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except: pass
    
    if not lig_coords: return "ERROR"

    # 3. Calculate Distances
    dists = distance.cdist(lig_coords, ghost_coords)
    min_dist = np.min(dists)
    
    # Threshold: 1.5 A (Hard Clash)
    return min_dist

def check_admet(ligand_path):
    """
    Calculates Lipinski's properties using RDKit.
    Note: PDBQT is not ideal for RDKit, usually need SDF/SMILES. 
    We will try to sanitize it best we can or use basic counts.
    """
    # Simply reading PDBQT with naive RDKit often fails. 
    # For robust pipeline, we assume the user might provide a better format, 
    # or we try to infer from the text block.
    # As a fallback, we skip if read fails.
    mol = Chem.MolFromPDBFile(str(ligand_path), removeHs=False)
    if not mol:
        # Try finding a matching SMILES or SDF if it helps? 
        return None

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    return mw, logp


def check_environment():
    """Verifies that all necessary files exist."""
    required = [VINA_EXEC, WT_RECEPTOR, MUT_RECEPTOR]
    missing = [str(f) for f in required if not f.exists()]
    if missing:
        print("Error: Missing required files:")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)
    
    # Make sure vina is executable
    if not os.access(VINA_EXEC, os.X_OK):
        print(f"Making {VINA_EXEC.name} executable...")
        try:
            os.chmod(VINA_EXEC, 0o755)
        except Exception as e:
            print(f"Warning: Could not set executable permissions: {e}")

def get_affinity(pdbqt_output):
    """Parses the best affinity from Vina output file."""
    best_affinity = 100.0
    found = False
    try:
        with open(pdbqt_output, 'r') as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT"):
                    parts = line.split()
                    if len(parts) >= 4:
                        affinity = float(parts[3])
                        if affinity < best_affinity:
                            best_affinity = affinity
                        found = True
    except FileNotFoundError:
        return None
    
    return best_affinity if found else None

def run_vina(ligand_path, receptor_path, output_path, seed=42):
    """Runs Vina docking."""
    cmd = [
        str(VINA_EXEC),
        "--receptor", str(receptor_path),
        "--ligand", str(ligand_path),
        "--center_x", str(BOX_CONFIG["center_x"]),
        "--center_y", str(BOX_CONFIG["center_y"]),
        "--center_z", str(BOX_CONFIG["center_z"]),
        "--size_x", str(BOX_CONFIG["size_x"]),
        "--size_y", str(BOX_CONFIG["size_y"]),
        "--size_z", str(BOX_CONFIG["size_z"]),
        "--exhaustiveness", str(BOX_CONFIG["exhaustiveness"]),
        "--num_modes", str(BOX_CONFIG["num_modes"]),
        "--scoring", "vinardo",
        "--seed", str(seed),
        "--out", str(output_path),
        "--cpu", "4" # Use multithreading if available
    ]
    
    # Run quietly
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print(f"Error: Docking failed for {receptor_path.name}")
        return False
    return True

def analyze_ligand(ligand_path):
    """Runs the full analysis pipeline for a single ligand."""
    ligand_file = Path(ligand_path)
    if not ligand_file.exists():
        print(f"Error: Ligand file {ligand_file} not found.")
        return

    print(f"\nAnalyzing Candidate: {ligand_file.name}")
    print("-" * 50)

    # Temporary output files
    out_wt = ligand_file.with_suffix('.wt_docked.pdbqt')
    out_mut = ligand_file.with_suffix('.mut_docked.pdbqt')

    # 1. Dock into Wild-Type
    print("  1. Docking into Wild-Type (9KNZ)... ", end="", flush=True)
    if run_vina(ligand_file, WT_RECEPTOR, out_wt):
        wt_score = get_affinity(out_wt)
        print(f"Done. (Affinity: {wt_score} kcal/mol)")
    else:
        print("Failed.")
        return

    # 2. Dock into Mutant (W730A)
    print("  2. Docking into Mutant (W730A)...   ", end="", flush=True)
    if run_vina(ligand_file, MUT_RECEPTOR, out_mut):
        mut_score = get_affinity(out_mut)
        print(f"Done. (Affinity: {mut_score} kcal/mol)")
    else:
        print("Failed.")
        return

    # 3. Calculate Delta
    delta = mut_score - wt_score
    print("-" * 50)
    print(f"  > Resistance Delta (Mutant - WT): {delta:+.2f} kcal/mol")
    
    # 3b. Ghost Clash Check (Rigorous Physics)
    clash_dist = check_ghost_clash(out_mut)
    if isinstance(clash_dist, float):
        print(f"  > Dist to Ghost Sidechain: {clash_dist:.2f} Å", end="")
        if clash_dist < 1.5:
            print(" [CLASH DETECTED]")
        else:
            print(" [SAFE]")

    # 4. Verdict
    print("-" * 50)
    
    passed_filters = True
    
    # Filter 1: Binding Strength
    if wt_score > -6.5:
        print("  VERDICT: REJECTED (Too weak to be a drug)")
        passed_filters = False
    
    # Filter 2: Resistance Resilience
    elif delta > 0.5:
        print("  VERDICT: REJECTED (Vulnerable to W730A mutation)")
        passed_filters = False
        
    # Filter 3: Ghost Clash
    elif isinstance(clash_dist, float) and clash_dist < 1.5:
        print("  VERDICT: REJECTED (Clashes with Wild-Type structure - Vacuum Hole Fallacy)")
        passed_filters = False

    if passed_filters:
        print("  VERDICT: ACCEPTED (Potent & Resilient Candidate!)")
        
        # Optional ADMET
        admet = check_admet(ligand_file) # Try original file first
        if not admet: admet = check_admet(out_wt) # Try pdbqt
        
        if admet:
            mw, logp = admet
            print(f"  ADMET Info: MW={mw:.1f}, LogP={logp:.1f}")
            if mw > 500: print("    Caution: Violates Lipinski MW rule")
            
        print("  Recommendation: Proceed to In-Vitro Validation.")
    
    # Cleanup
    if out_wt.exists(): out_wt.unlink()
    if out_mut.exists(): out_mut.unlink()

# =============================================================================
# MAIN
# =============================================================================
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Nipah Virus Inhibitor Discovery Pipeline")
    parser.add_argument("ligand", help="Path to the ligand PDBQT file to test")
    args = parser.parse_args()

    check_environment()
    analyze_ligand(args.ligand)
