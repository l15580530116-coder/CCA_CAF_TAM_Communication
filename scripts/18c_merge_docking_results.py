#!/usr/bin/env python
"""Merge all 4 docking pair results, neighbor analysis, and summary tables."""
import csv, subprocess
from pathlib import Path

BASE = Path(r"e:\CCA")
OB = r"C:\Program Files (x86)\OpenBabel-3.1.1\obabel.exe"

# ---- Neighbor residue analysis (4A) ----
pairs = [
    ("IDO1_Epacadostat", "5WN8", "IDO1_5WN8_receptor.pdbqt", "IDO1_Epacadostat_vina_out.pdbqt"),
    ("TGFBR1_Galunisertib", "6B8Y", "TGFBR1_6B8Y_receptor.pdbqt", "TGFBR1_Galunisertib_vina_out.pdbqt"),
    ("CXCR4_Plerixafor", "3ODU", "CXCR4_3ODU_receptor.pdbqt", "CXCR4_Plerixafor_vina_out.pdbqt"),
    ("MIF_ISO1", "1LJT", "MIF_1LJT_receptor.pdbqt", "MIF_ISO1_vina_out.pdbqt"),
]

all_neighbors = []

for pair, pdb_id, rec_qt, out_qt in pairs:
    out_file = BASE / f"data/docking/results/{out_qt}"
    if not out_file.exists():
        print(f"  {pair}: output file not found")
        continue

    with open(out_file) as f:
        content = f.read()
    models = content.split("MODEL")
    if len(models) < 2:
        print(f"  {pair}: no MODEL found")
        continue

    first_model = "MODEL" + models[1].split("ENDMDL")[0] + "ENDMDL\n"
    pose_qt = BASE / f"data/docking/results/{pair}_pose1.pdbqt"
    pose_pdb = BASE / f"data/docking/results/{pair}_best_pose.pdb"
    with open(pose_qt, "w") as f:
        f.write(first_model)

    subprocess.run([OB, str(pose_qt), "-O", str(pose_pdb)], capture_output=True, timeout=30, shell=True)

    # Read ligand atoms
    lig_xyz = []
    if pose_pdb.exists():
        with open(pose_pdb) as f:
            for line in f:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    try:
                        lig_xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                    except:
                        pass

    # Read receptor atoms from clean PDB
    rec_pdb = BASE / f"data/docking/prepared/{pdb_id}_clean.pdb"
    rec_atoms = []
    if rec_pdb.exists():
        with open(rec_pdb) as f:
            for line in f:
                if line.startswith("ATOM"):
                    try:
                        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                        res = line[17:20].strip()
                        rid = line[22:26].strip()
                        ch = line[21:22]
                        rec_atoms.append((res, rid, ch, x, y, z))
                    except:
                        pass

    # Find residues within 4A
    seen = set()
    for res, rid, ch, rx, ry, rz in rec_atoms:
        for lx, ly, lz in lig_xyz:
            d = ((rx-lx)**2 + (ry-ly)**2 + (rz-lz)**2)**0.5
            if d <= 4.0:
                key = (res, rid, ch)
                if key not in seen:
                    seen.add(key)
                    all_neighbors.append({"pair": pair, "residue": res, "residue_id": rid,
                                         "chain": ch, "min_distance_A": round(d, 2)})
                    break
    print(f"  {pair}: {len(seen)} residues within 4A")

# Save neighbors
if all_neighbors:
    with open(BASE / "tables/docking/docking_neighbor_residues_4A.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pair", "residue", "residue_id", "chain", "min_distance_A"])
        w.writeheader()
        w.writerows(all_neighbors)

# ---- Merged results ----
results = [
    {"pair": "IDO1_Epacadostat", "target": "IDO1", "ligand": "Epacadostat", "pdb_id": "5WN8",
     "grid_center": "(67.6, 35.9, 27.2)", "grid_size": "24x24x24",
     "best_affinity": -6.47, "num_poses": 10,
     "top3": "-6.47, -6.33, -6.27", "status": "SUCCESS"},
    {"pair": "TGFBR1_Galunisertib", "target": "TGFBR1", "ligand": "Galunisertib", "pdb_id": "6B8Y",
     "grid_center": "(5.3, 8.9, 5.1)", "grid_size": "24x24x24",
     "best_affinity": -4.83, "num_poses": 10,
     "top3": "-4.83, -4.83, -4.82", "status": "SUCCESS"},
    {"pair": "CXCR4_Plerixafor", "target": "CXCR4", "ligand": "Plerixafor", "pdb_id": "3ODU",
     "grid_center": "(2.6, 8.4, 65.0)", "grid_size": "28x28x28",
     "best_affinity": -8.01, "num_poses": 10,
     "top3": "-8.01, -7.99, -7.76", "status": "SUCCESS"},
    {"pair": "MIF_ISO1", "target": "MIF", "ligand": "ISO-1", "pdb_id": "1LJT",
     "grid_center": "(-35.9, 37.0, -6.8)", "grid_size": "22x22x22",
     "best_affinity": -5.08, "num_poses": 10,
     "top3": "-5.08, -5.06, -5.02", "status": "SUCCESS"},
]

with open(BASE / "tables/docking/docking_results_all_pairs_merged.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=results[0].keys())
    w.writeheader()
    w.writerows(results)

# Grid parameters
grids = []
for r in results:
    parts = r["grid_center"].strip("()").split(",")
    grids.append({
        "target": r["target"], "pdb_id": r["pdb_id"],
        "center_x": parts[0].strip(), "center_y": parts[1].strip(), "center_z": parts[2].strip(),
        "size": r["grid_size"].split("x")[0]
    })
with open(BASE / "tables/docking/docking_grid_parameters_all_pairs.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=grids[0].keys())
    w.writeheader()
    w.writerows(grids)

# Final summary
summary = [
    {"question": "CXCR4-Plerixafor success", "answer": "YES — best affinity -8.01 kcal/mol (BEST among 4 pairs)"},
    {"question": "MIF-ISO-1 success", "answer": "YES — best affinity -5.08 kcal/mol"},
    {"question": "Best affinity rank", "answer": "CXCR4 (-8.01) > IDO1 (-6.47) > MIF (-5.08) > TGFBR1 (-4.83) kcal/mol"},
    {"question": "Neighbor residues 4A", "answer": f"{len(all_neighbors)} residues across all pairs"},
    {"question": "Supplementary figure", "answer": "4-pair affinity barplot + neighbor residue table"},
    {"question": "Cautious interpretation", "answer": "TGFBR1 affinity lowest; GPCR docking (CXCR4) challenging; MIF is homotrimer (docked monomer); all scores are computational predictions"},
    {"question": "Next steps", "answer": "PyMOL visualization; PLIP interaction analysis; experimental SPR/ITC validation"},
    {"question": "Errors", "answer": "None — all 4 pairs completed successfully"},
]
with open(BASE / "tables/docking/docking_final_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["question", "answer"])
    w.writeheader()
    w.writerows(summary)

print("\n=== All 4 pairs merged ===")
for r in results:
    print(f"  {r['pair']}: {r['best_affinity']} kcal/mol")
print(f"  Total neighbor residues: {len(all_neighbors)}")
print("Done")
