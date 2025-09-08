#!/usr/bin/env python3
"""
Pipeline Smoke Test: Single Compound End-to-End Validation (v7 Protocol)

Runs ONE compound through the entire docking pipeline to verify
everything works before committing to the full 20-compound screen.

CRITICAL CHANGES from original:
- Uses Meeko (not OpenBabel) for ligand PDBQT
- exhaustiveness=16 (not 4)
- Deterministic: seed=42
- Includes RMSD redocking validation
- Bond perception for crystal ligand
"""

import os
import json
import subprocess
import sys
from pathlib import Path
import numpy as np

# Set up paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
LIGANDS_DIR = DATA_DIR / "ligands"
CONFIG_DIR = BASE_DIR / "config"
RESULTS_DIR = BASE_DIR / "results"
SCRIPTS_DIR = BASE_DIR / "scripts"
VINA_PATH = SCRIPTS_DIR / "vina"

# Ensure directories exist
DATA_DIR.mkdir(exist_ok=True)
LIGANDS_DIR.mkdir(exist_ok=True)
CONFIG_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

# Global settings (v7 protocol)
RANDOM_SEED = 42
EXHAUSTIVENESS = 16  # NOT 4! Critical for accuracy
BOX_SIZE = 22  # Angstroms

# Status tracking
status = {}

def print_step(step_num, description):
    print(f"\n{'='*60}")
    print(f"Step {step_num}: {description}")
    print('='*60)


# =============================================================================
# Step 1: Download and Clean PDB 9KNZ
# =============================================================================
print_step(1, "Download and Clean PDB 9KNZ")

import requests
from Bio.PDB import PDBParser, PDBIO, Select
import gemmi

# Download CIF file (PDB format not available for EM structures)
cif_url = "https://files.rcsb.org/download/9KNZ.cif"
cif_path = DATA_DIR / "9KNZ_raw.cif"
pdb_raw_path = DATA_DIR / "9KNZ_raw.pdb"

print(f"Downloading PDB 9KNZ (CIF format) from RCSB...")
response = requests.get(cif_url)
if response.status_code == 200:
    with open(cif_path, 'w') as f:
        f.write(response.text)
    print(f"  Downloaded CIF to: {cif_path}")

    # Convert CIF to PDB using gemmi
    print(f"  Converting CIF to PDB format...")
    doc = gemmi.cif.read(str(cif_path))
    block = doc.sole_block()
    structure = gemmi.make_structure_from_block(block)
    structure.write_pdb(str(pdb_raw_path))
    print(f"  Converted to: {pdb_raw_path}")
else:
    print(f"  ERROR: Failed to download CIF (status {response.status_code})")
    sys.exit(1)

# Clean: isolate Chain A, remove waters/ions
class ChainASelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"
    def accept_residue(self, residue):
        return residue.id[0] == " "  # Exclude heteroatoms/waters

parser = PDBParser(QUIET=True)
structure = parser.get_structure("9KNZ", str(pdb_raw_path))

# Check what chains are available
chains = [c.id for model in structure for c in model]
print(f"  Available chains: {chains}")

# If Chain A doesn't exist, use first available chain
if "A" not in chains:
    first_chain = chains[0] if chains else None
    print(f"  WARNING: Chain A not found, using Chain {first_chain}")
    class ChainASelect(Select):
        def accept_chain(self, chain):
            return chain.id == first_chain
        def accept_residue(self, residue):
            return residue.id[0] == " "

io = PDBIO()
io.set_structure(structure)
clean_pdb_path = DATA_DIR / "9KNZ_clean.pdb"
io.save(str(clean_pdb_path), ChainASelect())

# Verify
line_count = sum(1 for _ in open(clean_pdb_path))
print(f"  Saved cleaned structure: {clean_pdb_path}")
print(f"  Line count: {line_count}")
status["step1_structure"] = line_count > 1000
print(f"  CHECK: {'PASS' if status['step1_structure'] else 'FAIL'} (expected >1000 lines)")


# =============================================================================
# Step 2: Extract ERDRP-0519 Ligand (for redocking validation)
# =============================================================================
print_step(2, "Extract ERDRP-0519 Ligand")

class LigandSelect(Select):
    def accept_residue(self, residue):
        # Include heteroatoms but exclude common ions and water
        return residue.id[0] != " " and residue.resname not in ["HOH", "NA", "CL", "MG", "ZN", "CA", "K"]

# Re-parse the raw structure to get ligands
structure = parser.get_structure("9KNZ", str(pdb_raw_path))
io = PDBIO()
io.set_structure(structure)
ligand_path = DATA_DIR / "ERDRP_original.pdb"
io.save(str(ligand_path), LigandSelect())

# Verify and show what ligands we found
line_count = sum(1 for _ in open(ligand_path))
print(f"  Extracted ligand to: {ligand_path}")
print(f"  Line count: {line_count}")

# List unique ligand residue names
ligand_resnames = set()
with open(ligand_path) as f:
    for line in f:
        if line.startswith(("HETATM", "ATOM")):
            resname = line[17:20].strip()
            ligand_resnames.add(resname)
print(f"  Ligand residues found: {ligand_resnames}")

status["step2_ligand"] = line_count > 10
print(f"  CHECK: {'PASS' if status['step2_ligand'] else 'FAIL'} (expected >10 lines)")


# =============================================================================
# Step 3: Calculate Docking Box from Crystal Ligand
# =============================================================================
print_step(3, "Calculate Docking Box from Crystal Ligand")

coords = []
ref_coords = []  # For RMSD calculation later

with open(ligand_path) as f:
    for line in f:
        if line.startswith(("HETATM", "ATOM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
                # Store heavy atom coords for RMSD (skip H)
                elem = line[76:78].strip() if len(line) > 76 else line[12:14].strip()
                if elem not in ["H", "D"]:
                    ref_coords.append([x, y, z])
            except ValueError:
                pass

if coords:
    center_x = sum(c[0] for c in coords) / len(coords)
    center_y = sum(c[1] for c in coords) / len(coords)
    center_z = sum(c[2] for c in coords) / len(coords)
    print(f"  Ligand center: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f})")
    
    # Save reference coords for RMSD
    np.save(DATA_DIR / "ERDRP_ref_coords.npy", np.array(ref_coords))
    print(f"  Saved {len(ref_coords)} heavy atom reference coordinates")
else:
    print("  ERROR: No ligand coordinates found!")
    sys.exit(1)

docking_config = {
    "center_x": round(center_x, 3),
    "center_y": round(center_y, 3),
    "center_z": round(center_z, 3),
    "size_x": BOX_SIZE,
    "size_y": BOX_SIZE,
    "size_z": BOX_SIZE,
    "exhaustiveness": EXHAUSTIVENESS,
    "seed": RANDOM_SEED
}

config_path = CONFIG_DIR / "docking_box.json"
with open(config_path, "w") as f:
    json.dump(docking_config, f, indent=2)

print(f"  Saved config: {config_path}")
print(f"  Config: {docking_config}")

status["step3_box"] = config_path.exists()
print(f"  CHECK: {'PASS' if status['step3_box'] else 'FAIL'}")


# =============================================================================
# Step 4: Prepare Receptor PDBQT Files (OpenBabel OK for receptor)
# =============================================================================
print_step(4, "Prepare Receptor PDBQT Files")

wt_pdbqt = DATA_DIR / "9KNZ_clean.pdbqt"

# Add hydrogens at pH 7.4, then convert to PDBQT
print("  Adding hydrogens at pH 7.4...")
wt_h_pdb = DATA_DIR / "9KNZ_clean_H.pdb"
cmd1 = ["obabel", str(clean_pdb_path), "-O", str(wt_h_pdb), "-p", "7.4"]
result1 = subprocess.run(cmd1, capture_output=True, text=True)

print("  Converting to PDBQT...")
cmd2 = ["obabel", str(wt_h_pdb), "-O", str(wt_pdbqt), "-xr"]
result2 = subprocess.run(cmd2, capture_output=True, text=True)

if wt_pdbqt.exists():
    print(f"  Created: {wt_pdbqt}")
else:
    print(f"  ERROR: Failed to create receptor PDBQT")

status["step4_receptor"] = wt_pdbqt.exists()
print(f"  CHECK: {'PASS' if status['step4_receptor'] else 'FAIL'}")


# =============================================================================
# Step 5: Prepare Crystal Ligand with Bond Perception (CRITICAL)
# =============================================================================
print_step(5, "Prepare Crystal Ligand (Bond Perception + Meeko)")

# from openbabel import openbabel as ob # Removed to avoid build errors
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

print("  Step 5a: Fixing PDB format...")
fixed_lines = []
atom_num = 0

with open(ligand_path, "r") as f:
    for line in f:
        if line.startswith("HETATM") or line.startswith("ATOM"):
            atom_num += 1
            parts = line.split()
            if len(parts) >= 8:
                name = parts[2]
                elem = parts[-1] if len(parts[-1]) <= 2 else name[0]
                x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
                fixed = f"HETATM{atom_num:>5} {name:<4} LIG A   1    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {elem:>2}\n"
                fixed_lines.append(fixed)

fixed_lines.append("END\n")
fixed_pdb = DATA_DIR / "ERDRP_fixed.pdb"
with open(fixed_pdb, "w") as f:
    f.writelines(fixed_lines)
print(f"  Fixed PDB: {atom_num} atoms")

print("  Step 5b: Bond perception (using obabel CLI)...")
sdf_path = DATA_DIR / "ERDRP_with_bonds.sdf"

# Use obabel CLI to convert PDB -> SDF (automatically handles bond perception and Hydrogens)
# -i pdb: input format
# -o sdf: output format
# -h: add hydrogens (pH 7.4 is default for protein, but -h works for ligands too)
cmd = ["obabel", "-ipdb", str(fixed_pdb), "-osdf", "-O", str(sdf_path), "-h"]
result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode != 0:
    print(f"  ERROR: obabel CLI failed: {result.stderr}")
    sys.exit(1)
    
print(f"    Converted to SDF: {sdf_path}")

print("  Step 5c: RDKit validation...")
rdmol = Chem.MolFromMolFile(str(sdf_path), removeHs=False)
if rdmol is None:
    print("  ERROR: RDKit failed to read SDF!")
    sys.exit(1)

n_frags = len(Chem.GetMolFrags(rdmol))
print(f"    Fragments: {n_frags}")
if n_frags != 1:
    print(f"  ERROR: Bond perception FAILED! Molecule has {n_frags} fragments!")
    sys.exit(1)

print("  Step 5d: Meeko PDBQT preparation...")
preparator = MoleculePreparation()
mol_setups = preparator.prepare(rdmol)
pdbqt_string = PDBQTWriterLegacy.write_string(mol_setups[0])[0]

erdrp_pdbqt = LIGANDS_DIR / "ERDRP.pdbqt"
with open(erdrp_pdbqt, "w") as f:
    f.write(pdbqt_string)
print(f"  Created: {erdrp_pdbqt}")


status["step5_crystal_ligand"] = erdrp_pdbqt.exists() and n_frags == 1
print(f"  CHECK: {'PASS' if status['step5_crystal_ligand'] else 'FAIL'}")


# =============================================================================
# Step 6: Redocking Validation (RMSD < 2.0 Å)
# =============================================================================
print_step(6, "Redocking Validation (RMSD must be < 2.0 Å)")

with open(config_path) as f:
    config = json.load(f)

redock_output = RESULTS_DIR / "ERDRP_redocked.pdbqt"

cmd = [
    str(VINA_PATH),
    "--receptor", str(wt_pdbqt),
    "--ligand", str(erdrp_pdbqt),
    "--center_x", str(config["center_x"]),
    "--center_y", str(config["center_y"]),
    "--center_z", str(config["center_z"]),
    "--size_x", str(config["size_x"]),
    "--size_y", str(config["size_y"]),
    "--size_z", str(config["size_z"]),
    "--exhaustiveness", str(EXHAUSTIVENESS),
    "--seed", str(RANDOM_SEED),
    "--out", str(redock_output)
]

print(f"  Running Vina redocking (exhaustiveness={EXHAUSTIVENESS}, seed={RANDOM_SEED})...")
result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)

# Parse RMSD from docked output
def parse_pdbqt_heavy_atoms(filename, model=1):
    """Extract heavy atom coordinates from PDBQT."""
    coords = []
    current_model = 0
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = int(line.split()[1])
            if line.startswith(("ATOM", "HETATM")):
                if current_model == model or current_model == 0:
                    atom_type = line[77:79].strip() if len(line) > 77 else ""
                    if atom_type not in ["H", "HD", "HS"]:
                        try:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            coords.append([x, y, z])
                        except:
                            pass
            if line.startswith("ENDMDL") and current_model == model:
                break
    return np.array(coords)

def calculate_rmsd(coords1, coords2):
    """Calculate RMSD between two coordinate sets."""
    min_len = min(len(coords1), len(coords2))
    if min_len == 0:
        return 999.0
    diff = coords1[:min_len] - coords2[:min_len]
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

# Calculate RMSD
ref_coords = np.load(DATA_DIR / "ERDRP_ref_coords.npy")
docked_coords = parse_pdbqt_heavy_atoms(redock_output, model=1)
rmsd = calculate_rmsd(ref_coords, docked_coords)

print(f"  Reference atoms: {len(ref_coords)}")
print(f"  Docked atoms: {len(docked_coords)}")
print(f"  RMSD: {rmsd:.3f} Å")

if rmsd <= 2.0:
    print(f"  VALIDATION PASSED: RMSD ≤ 2.0 Å")
    status["step6_rmsd"] = True
else:
    print(f"  WARNING: RMSD > 2.0 Å - docking setup may have issues!")
    status["step6_rmsd"] = False

print(f"  CHECK: {'PASS' if status['step6_rmsd'] else 'FAIL'}")


# =============================================================================
# Step 7: Create W730A Mutant Structure
# =============================================================================
print_step(7, "Create W730A Mutant Structure")

TARGET_TRP = 730  # Nipah equivalent of Measles W671

structure = parser.get_structure("WT", str(clean_pdb_path))
trp_residues = []
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.resname == "TRP":
                trp_residues.append((chain.id, residue.id[1], residue.resname))

print(f"  TRP residues found: {[r[1] for r in trp_residues]}")
print(f"  Target for mutation: W{TARGET_TRP}")

# Verify residue 730 is Tryptophan
verified = False
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[1] == TARGET_TRP:
                if residue.resname == "TRP":
                    print(f"  VERIFIED: Residue {TARGET_TRP} is TRP")
                    verified = True
                else:
                    print(f"  ERROR: Residue {TARGET_TRP} is {residue.resname}, not TRP!")

if not verified:
    print(f"  WARNING: Residue {TARGET_TRP} not found as TRP!")
    # Try to find closest TRP
    if trp_residues:
        closest = min(trp_residues, key=lambda x: abs(x[1] - TARGET_TRP))
        print(f"  Using closest TRP: W{closest[1]}")
        TARGET_TRP = closest[1]
        verified = True

# Perform mutation: W -> A
ala_atoms = ["N", "CA", "C", "O", "CB"]
mutated = False

structure = parser.get_structure("WT", str(clean_pdb_path))
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[1] == TARGET_TRP and residue.resname == "TRP":
                atoms_to_remove = [atom.id for atom in residue if atom.id not in ala_atoms]
                for atom_id in atoms_to_remove:
                    residue.detach_child(atom_id)
                residue.resname = "ALA"
                mutated = True
                print(f"  Mutated W{TARGET_TRP} -> A{TARGET_TRP}")

io = PDBIO()
io.set_structure(structure)
mutant_path = DATA_DIR / "9KNZ_W730A.pdb"
io.save(str(mutant_path))

# Prepare mutant PDBQT
mut_h_pdb = DATA_DIR / "9KNZ_W730A_H.pdb"
mut_pdbqt = DATA_DIR / "9KNZ_W730A.pdbqt"
subprocess.run(["obabel", str(mutant_path), "-O", str(mut_h_pdb), "-p", "7.4"], capture_output=True)
subprocess.run(["obabel", str(mut_h_pdb), "-O", str(mut_pdbqt), "-xr"], capture_output=True)

print(f"  Saved mutant: {mutant_path}")
print(f"  Saved mutant PDBQT: {mut_pdbqt}")

status["step7_mutant"] = mutant_path.exists() and mut_pdbqt.exists() and mutated
print(f"  CHECK: {'PASS' if status['step7_mutant'] else 'FAIL'}")


# =============================================================================
# Step 8: Prepare Test Ligand (Ribavirin) with Meeko
# =============================================================================
print_step(8, "Prepare Test Ligand (Ribavirin) with Meeko")

from rdkit.Chem import AllChem

# Ribavirin SMILES
test_smiles = "C1=NC(=NN1C2C(C(C(O2)CO)O)O)C(=O)N"
test_name = "Ribavirin"

print(f"  Generating 3D structure for {test_name}...")
mol = Chem.MolFromSmiles(test_smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=RANDOM_SEED)
AllChem.MMFFOptimizeMolecule(mol)

# Use Meeko for PDBQT (NOT OpenBabel!)
print(f"  Converting to PDBQT with Meeko...")
preparator = MoleculePreparation()
mol_setups = preparator.prepare(mol)
pdbqt_string = PDBQTWriterLegacy.write_string(mol_setups[0])[0]

pdbqt_path = LIGANDS_DIR / f"{test_name}.pdbqt"
with open(pdbqt_path, "w") as f:
    f.write(pdbqt_string)
print(f"  Created: {pdbqt_path}")

status["step8_test_ligand"] = pdbqt_path.exists()
print(f"  CHECK: {'PASS' if status['step8_test_ligand'] else 'FAIL'}")


# =============================================================================
# Step 9: Dock Ribavirin to Wildtype
# =============================================================================
print_step(9, "Dock Ribavirin to Wildtype")

wt_output = RESULTS_DIR / "wt_docked.pdbqt"

cmd = [
    str(VINA_PATH),
    "--receptor", str(wt_pdbqt),
    "--ligand", str(pdbqt_path),
    "--center_x", str(config["center_x"]),
    "--center_y", str(config["center_y"]),
    "--center_z", str(config["center_z"]),
    "--size_x", str(config["size_x"]),
    "--size_y", str(config["size_y"]),
    "--size_z", str(config["size_z"]),
    "--exhaustiveness", str(EXHAUSTIVENESS),
    "--seed", str(RANDOM_SEED),
    "--out", str(wt_output)
]

print(f"  Running Vina (exhaustiveness={EXHAUSTIVENESS}, seed={RANDOM_SEED})...")
result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)

# Parse affinity
result_wt = None
for line in result.stdout.split('\n'):
    if line.strip().startswith('1'):
        parts = line.split()
        if len(parts) >= 2:
            try:
                result_wt = float(parts[1])
                break
            except ValueError:
                pass

print(f"  Wildtype affinity: {result_wt} kcal/mol")
status["step9_wt_dock"] = result_wt is not None
print(f"  CHECK: {'PASS' if status['step9_wt_dock'] else 'FAIL'}")


# =============================================================================
# Step 10: Dock Ribavirin to Mutant
# =============================================================================
print_step(10, "Dock Ribavirin to Mutant (W730A)")

mut_output = RESULTS_DIR / "mut_docked.pdbqt"

cmd = [
    str(VINA_PATH),
    "--receptor", str(mut_pdbqt),
    "--ligand", str(pdbqt_path),
    "--center_x", str(config["center_x"]),
    "--center_y", str(config["center_y"]),
    "--center_z", str(config["center_z"]),
    "--size_x", str(config["size_x"]),
    "--size_y", str(config["size_y"]),
    "--size_z", str(config["size_z"]),
    "--exhaustiveness", str(EXHAUSTIVENESS),
    "--seed", str(RANDOM_SEED),
    "--out", str(mut_output)
]

print(f"  Running Vina...")
result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)

result_mut = None
for line in result.stdout.split('\n'):
    if line.strip().startswith('1'):
        parts = line.split()
        if len(parts) >= 2:
            try:
                result_mut = float(parts[1])
                break
            except ValueError:
                pass

print(f"  Mutant affinity: {result_mut} kcal/mol")
status["step10_mut_dock"] = result_mut is not None

# Step 11: Results
print_step(11, "Calculate Delta Affinity")

import pandas as pd

if result_wt is not None and result_mut is not None:
    delta_affinity = result_mut - result_wt
    print(f"  WT: {result_wt:.2f}, Mut: {result_mut:.2f}, Delta: {delta_affinity:.2f} kcal/mol")
    
    results = pd.DataFrame([{"compound": "Ribavirin", "wt": result_wt, "mut": result_mut, "delta": delta_affinity}])
    results.to_csv(RESULTS_DIR / "smoke_test_results.csv", index=False)
    status["step11_results"] = True

# Summary
print("\n" + "="*60)
print("SMOKE TEST SUMMARY")
print("="*60)
for k, v in status.items():
    print(f"  {'PASS' if v else 'FAIL'}: {k}")
print("="*60)
