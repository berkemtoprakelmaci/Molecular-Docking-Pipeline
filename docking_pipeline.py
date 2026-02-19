# Copyright (C) 2026 Berkem Toprak Elmacı
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.

import os
import sys
import subprocess
import numpy as np

# SETTINGS
PDB_ID       = "1stp"          # PDB code of the protein
CHAIN        = "A"             # Chain for use
LIGAND_RESN  = "BTN"           # Ligand residue name (3-letter code)
PADDING      = 10.0            # Grid box extra space on each side (Angstrom)
EXHAUSTIVENESS = 8             # Vina exhaustiveness (higher, better but slower)
NUM_MODES    = 9               # How many poses to calculate

# Vina path (vina.exe for Windows, vina for Linux/Mac)
VINA_EXE = "vina.exe" if sys.platform == "win32" else "vina"

def run(cmd, description=""):
    """Run a command, stop if there is an error."""
    if description:
        print(f"  >> {description}")
    print(f"     Command: {cmd}")
    ret = subprocess.run(cmd, shell=True)
    if ret.returncode != 0:
        print(f"\n[ERROR] The following command failed (code {ret.returncode}):")
        print(f"  {cmd}")
        print("Pipeline stopped")
        sys.exit(1)
    print("     OK\n")


def step1_prepare_protein():
    """Download protein with PyMOL, remove water, add hydrogens, save."""
    print("=" * 55)
    print("STEP 1: Preparing protein and ligand (PyMOL)")
    print("=" * 55)

    pml_content = f"""
# Download protein
fetch {PDB_ID}, async=0

# Remove water first, then add hydrogens - correct order
remove resn HOH

# Add hydrogens
h_add

# Save protein (selected chain, polymer only)
save receptor.pdb, chain {CHAIN} and polymer

# Save ligand separately (coordinates needed for grid calculation)
save ligand_ref.pdb, resn {LIGAND_RESN}

quit
"""

    with open("prepare.pml", "w") as f:
        f.write(pml_content)

    run("pymol -c prepare.pml", "Preparing protein with PyMOL")

    # File check
    for file in ["receptor.pdb", "ligand_ref.pdb"]:
        if not os.path.exists(file) or os.path.getsize(file) == 0:
            print(f"[ERROR] {file} could not be created or is empty!")
            print("  - Is the PDB ID correct? Is the ligand residue name correct?")
            sys.exit(1)

    print("  Protein and ligand PDB files are ready.\n")


def step2_calculate_grid():
    """
    Read atom coordinates from ligand reference file,
    automatically calculate grid center and size.
    Returns: (cx, cy, cz, sx, sy, sz)
    """
    print("=" * 55)
    print("STEP 2: Calculating automatic grid (ligand coordinates)")
    print("=" * 55)

    coordinates = []

    with open("ligand_ref.pdb", "r") as f:
        for line in f:
            # Read ATOM or HETATM lines
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coordinates.append([x, y, z])
                except ValueError:
                    continue

    if not coordinates:
        print("[ERROR] No coordinates found in ligand_ref.pdb!")
        print(f"  Is the ligand residue name correct? ({LIGAND_RESN})")
        sys.exit(1)

    coords = np.array(coordinates)

    # Center: average of coordinates
    center = coords.mean(axis=0)

    # Size: min-max difference + padding (on both sides)
    size = (coords.max(axis=0) - coords.min(axis=0)) + 2 * PADDING

    # Minimum 15 Angstrom
    size = np.maximum(size, 15.0)

    cx, cy, cz = center
    sx, sy, sz = size

    print(f"  Atom count     : {len(coordinates)}")
    print(f"  Grid center    : X={cx:.2f}  Y={cy:.2f}  Z={cz:.2f}")
    print(f"  Grid size      : X={sx:.2f}  Y={sy:.2f}  Z={sz:.2f}")
    print()

    return cx, cy, cz, sx, sy, sz


def step3_convert_formats():
    """Convert PDB -> PDBQT with Open Babel."""
    print("=" * 55)
    print("STEP 3: Converting formats (Open Babel)")
    print("=" * 55)

    # Receptor: rigid, obabel handles Gasteiger charge assignment
    # -xr = rigid (lock all rotatable bonds)
    # -h  = add implicit hydrogens (no problem if already added)
    
    run(
        "obabel receptor.pdb -O receptor.pdbqt -xr",
        "Receptor PDB -> PDBQT"
    )

    run(
        "obabel ligand_ref.pdb -O ligand.pdbqt -h -p 7.4",
        "Ligand PDB -> PDBQT (pH 7.4)"
    )

    for file in ["receptor.pdbqt", "ligand.pdbqt"]:
        if not os.path.exists(file) or os.path.getsize(file) == 0:
            print(f"[ERROR] {file} could not be created! Is Open Babel installed?")
            sys.exit(1)

    print("  PDBQT files are ready.\n")


def step4_run_docking(cx, cy, cz, sx, sy, sz):
    """Run docking with Vina."""
    print("=" * 55)
    print("STEP 4: Running docking (AutoDock Vina)")
    print("=" * 55)

    vina_command = (
        f"{VINA_EXE} "
        f"--receptor receptor.pdbqt "
        f"--ligand ligand.pdbqt "
        f"--center_x {cx:.3f} --center_y {cy:.3f} --center_z {cz:.3f} "
        f"--size_x {sx:.3f} --size_y {sy:.3f} --size_z {sz:.3f} "
        f"--exhaustiveness {EXHAUSTIVENESS} "
        f"--num_modes {NUM_MODES} "
        f"--out result.pdbqt"
    )

    # --log was removed in Vina 1.2.7
    # capturing output with Python and saving to file
    
    print("  >> Starting Vina docking")
    print(f"     Command: {vina_command}")
    ret = subprocess.run(vina_command, shell=True, capture_output=True, text=True)
    output = ret.stdout + ret.stderr
    print(output)
    with open("result_log.txt", "w", encoding="utf-8") as f:
        f.write(output)
    if ret.returncode != 0:
        print(f"\n[ERROR] Vina failed (code {ret.returncode})")
        sys.exit(1)
    print("     OK\n")
    print("  Docking completed!\n")


def step5_show_results():
    """Read and print results from log file."""
    print("=" * 55)
    print("STEP 5: Results")
    print("=" * 55)

    if not os.path.exists("result_log.txt"):
        print("  Log file not found, check result.pdbqt file.")
        return

    with open("result_log.txt", "r") as f:
        content = f.read()

    print(content)

    print("\nCreated files:")
    for file in ["receptor.pdbqt", "ligand.pdbqt", "result.pdbqt", "result_log.txt"]:
        if os.path.exists(file):
            size = os.path.getsize(file)
            print(f"  {file:25s}  ({size} bytes)")

# MAIN PROGRAM

if __name__ == "__main__":

    step1_prepare_protein()
    cx, cy, cz, sx, sy, sz = step2_calculate_grid()
    step3_convert_formats()
    step4_run_docking(cx, cy, cz, sx, sy, sz)
    step5_show_results()
    print("Next steps:")
    print("  1. Check binding affinity scores in result_log.txt")
    print("  2. Open result.pdbqt in PyMOL to visualize poses")
    print("  3. Load the best pose together with receptor.pdb and analyze")
