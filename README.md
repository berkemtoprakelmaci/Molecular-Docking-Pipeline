# Molecular Docking Pipeline

A simple Python pipeline that automatically downloads a protein, prepares it, calculates the grid, and runs docking.

## Requirements
- [Miniconda3] 
- [PyMOL](https://pymol.org) (command line version)
- [Open Babel](https://openbabel.org)
- [AutoDock Vina 1.2+](https://github.com/ccsb-scripps/AutoDock-Vina/releases)
- Python 3.10+

## Installation

```bash
conda create -n docking_project python=3.9 -y
conda activate docking_project
conda install -c conda-forge pymol-open-source openbabel -y
```
> AutoDock Vina must be installed separately and added to directory as `vina.exe` that consist `docking_pipeline.py`.
---

## Usage

1. Open `docking_pipeline.py` and edit the settings at the top:

```python
PDB_ID       = "1stp"          # PDB code of the protein
CHAIN        = "A"             # Chain for use
LIGAND_RESN  = "BTN"           # Ligand residue name (3-letter code)
PADDING      = 10.0            # Grid box extra space on each side (Angstrom)
EXHAUSTIVENESS = 8             # Vina exhaustiveness (higher, better but slower)
NUM_MODES    = 9               # How many poses to calculate
```

2. Run:
Open anaconda bash at the directory that contains `docking_pipeline.py` and `vina.exe`
```anaconda bash
conda activate docking_project
docking_pipeline.py
```

## Output Files

| File | Description |
| `receptor.pdbqt` | Prepared receptor |
| `ligand.pdbqt` | Prepared ligand |
| `result.pdbqt` | Docking poses (9 modes) |
| `result_log.txt` | Binding affinity scores |

## Visualizing Results

Open in PyMOL:
```
load receptor.pdb
load result.pdbqt
```
or
```
pymol receptor.pdb result.pdbqt
```
## Example Result (1stp / Biotin)

```
mode |   affinity (kcal/mol)
-----+-----------------------
   1 |   -7.109
   2 |   -6.153
   3 |   -6.007
```

> A score below -7 kcal/mol is generally considered good binding.


