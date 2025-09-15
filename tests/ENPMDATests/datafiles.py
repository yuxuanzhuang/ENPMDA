from pathlib import Path

# Resolve paths relative to this file's directory
_BASE = Path(__file__).resolve().parent

ensemble_ala_bonded = []
ensemble_ala_traj = []
ensemble_ala_top = []

for rep in range(1, 9):
    ensemble_ala_bonded.append(
        str(_BASE / f"datafile/ensemble_AlaDipeptide/rep{rep}/md.tpr")
    )
    ensemble_ala_traj.append(
        str(_BASE / f"datafile/ensemble_AlaDipeptide/rep{rep}/md.xtc")
    )
    ensemble_ala_top.append(
        str(_BASE / f"datafile/ensemble_AlaDipeptide/rep{rep}/start.pdb")
    )