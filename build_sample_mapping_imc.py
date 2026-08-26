"""
Build simulation_id -> sample_id -> patient_id mapping for a dataset whose PCMM driver
loops over a CSV in row order (see run_imc_wellmixed.jl / run_imc_spatial.jl: both loop
over `eachrow(df)` and push one monad per row, in file order, before any jobs execute --
same pattern already verified for htan_wellmixed).

Verified against the driver's own execution log ("Running simulation: N..." lines plus
"WARNING: Simulation N failed" lines) -- never inferred from filenames or folder naming,
per the lesson learned on htan_wellmixed's disproven <cell_positions><filename> approach.

Usage:
    python3 build_sample_mapping_imc.py <input_csv> <sample_id_column> <retry_map_json> <out_csv>

<retry_map_json> maps a failed original CSV-row position (1-based) to the simulation_id
that eventually succeeded for it, as a JSON object string, e.g.:
    '{"5": 49, "21": 50, "26": 51, "34": 52}'
Positions not listed are assumed to have succeeded on the first attempt (simulation_id
== row position). Pass '{}' if nothing ever failed for this dataset.

Patient ID is the leading letters+digits prefix (JHH368 from JHH368ROI1, HT056 from
HT056P1_S1PA, etc.) -- adjust PATIENT_REGEX if a future dataset uses something else.
"""
import sys
import json
import re
import pandas as pd
from pathlib import Path

PATIENT_REGEX = re.compile(r"^([A-Za-z]+\d+)")

if len(sys.argv) != 5:
    print("Usage: python3 build_sample_mapping_imc.py <input_csv> <sample_id_column> "
          "<retry_map_json> <out_csv>")
    sys.exit(1)

input_csv, sample_col, retry_json, out_csv = sys.argv[1:5]
retry_map = {int(k): v for k, v in json.loads(retry_json).items()}

df = pd.read_csv(input_csv)
if sample_col not in df.columns:
    print(f"ERROR: column '{sample_col}' not found in {input_csv}. "
          f"Available columns: {list(df.columns)}")
    sys.exit(1)

rows = []
for i, sample_id in enumerate(df[sample_col], start=1):
    simulation_id = retry_map.get(i, i)
    m = PATIENT_REGEX.match(str(sample_id))
    patient_id = m.group(1) if m else str(sample_id)
    rows.append({
        "simulation": str(simulation_id),
        "sample_id": sample_id,
        "patient_id": patient_id,
        "csv_row_position": i,
        "was_retried": i in retry_map,
    })

mapping = pd.DataFrame(rows).sort_values("simulation", key=lambda s: s.astype(int))
mapping.to_csv(out_csv, index=False)
print(f"Wrote {out_csv}: {len(mapping)} simulations mapped "
      f"({mapping['was_retried'].sum()} were retries)")
print(mapping.to_string(index=False))
