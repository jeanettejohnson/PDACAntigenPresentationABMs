"""
Build the simulation_id -> sample_id -> patient_id mapping for htan_wellmixed,
verified directly against PCMM's own execution logs (not filename inference):

  - slurm/logs/pcmm-driver-htan_wellmixed_19641529.out (Trial 1):
      queued 73 CSV rows in file order; "Max parallel sims: 1" means execution
      order == simulation_id order, so position N in the CSV == simulation N.
      Simulations 8, 9, 22, 43 FAILED (see WARNING lines in that log).

  - slurm/logs/pcmm-driver-htan_wellmixed_19647329.out (Trial 2):
      re-queued the same 73 rows; 69 matched already-successful sims and were
      skipped, the 4 that previously failed were re-run as brand-new
      simulation IDs 74, 75, 76, 77 (in that same relative order).

  Independently confirmed by the ic_cell_variations/ folder's csv provenance
  suffixes: ic_cell_variation_8_s74.csv, _9_s75.csv, _22_s76.csv, _43_s77.csv
  -- "_s<N>" there is exactly the retry's real simulation_id.

Usage (run from the PDACAntigenPresentationABMs repo root, where
assignmentsummary_HTAN_singlecell.csv lives):
    python3 build_sample_mapping_v2.py assignmentsummary_HTAN_singlecell.csv simulation_to_sample_id.csv
"""
import sys
import pandas as pd
from pathlib import Path

if len(sys.argv) != 3:
    print("Usage: python3 build_sample_mapping_v2.py <assignment_csv> <out_csv>")
    sys.exit(1)

assignment_csv = Path(sys.argv[1])
out_csv = Path(sys.argv[2])

# position (1-based, matching CSV row order) -> real simulation_id, for the
# four rows whose first attempt failed and got re-run under a new ID.
FAILED_RETRY_MAP = {8: 74, 9: 75, 22: 76, 43: 77}

df = pd.read_csv(assignment_csv)
if len(df) != 73:
    print(f"WARNING: expected 73 rows (matching the original Queuing log), got {len(df)}. "
          f"Double check this is the same assignment CSV used for the original run "
          f"before trusting the mapping below.")

rows = []
for i, sample_id in enumerate(df["sample_id"], start=1):
    simulation_id = FAILED_RETRY_MAP.get(i, i)
    patient_id = sample_id.split("P1_")[0] if "P1_" in sample_id else sample_id
    rows.append({
        "simulation": str(simulation_id),
        "sample_id": sample_id,
        "patient_id": patient_id,
        "csv_row_position": i,
        "was_retried": i in FAILED_RETRY_MAP,
    })

mapping = pd.DataFrame(rows).sort_values("simulation", key=lambda s: s.astype(int))
mapping.to_csv(out_csv, index=False)
print(f"Wrote {out_csv}: {len(mapping)} simulations mapped "
      f"({mapping['was_retried'].sum()} were retries)")
print(mapping.to_string(index=False))
