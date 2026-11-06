"""Generate instrument-neutral DC resistor-network acquisition plans.

No hardware commands are sent. Physical units are converted to the normalized
operator used by intervention.py. Recorded decimal readings are interval centers,
not exact physical observations. A real certificate also needs calibrated current,
shunt, relay/leakage and voltmeter uncertainty bounds.
"""
from __future__ import annotations

import argparse
import csv
import json
from fractions import Fraction as F
from pathlib import Path

from intervention import acquisition_plan


def plan(n, scale=F(1, 1000), current=F(1, 10000)):
    scale, current = F(scale), F(current)
    if scale <= 0 or current <= 0:
        raise ValueError("positive current and conductance scale required")
    rows = []
    for config in acquisition_plan(n):
        entries = sorted(config["loads"].items(), key=lambda x: int(x[0]))
        row = {"configuration": config["id"], "source_pad": "node_1",
               "sense_pad": "node_1", "current_A": str(current),
               "conductance_scale_S": str(scale), "readout_V": "",
               "voltage_absolute_error_V": "", "current_relative_error": "",
               "shunt_A_pad": "", "shunt_A_resistance_ohm": "",
               "shunt_A_absolute_error_S": "", "shunt_B_pad": "",
               "shunt_B_resistance_ohm": "", "shunt_B_absolute_error_S": "",
               "settled": "", "repeat_id": "", "temperature_C": ""}
        for label, (j, t) in zip(("A", "B"), entries):
            row[f"shunt_{label}_pad"] = f"node_{int(j)+1}"
            row[f"shunt_{label}_resistance_ohm"] = str(1/(F(t)*scale))
        rows.append(row)
    return rows


def voltage_interval(row):
    """Conservative normalization of voltage/current uncertainty only."""
    v = F(row["readout_V"])
    ev = F(row["voltage_absolute_error_V"])
    i = F(row["current_A"])
    ei = F(row["current_relative_error"])
    s = F(row["conductance_scale_S"])
    if ev < 0 or not 0 <= ei < 1 or i <= 0 or s <= 0 or v-ev <= 0:
        raise ValueError("positive signal and valid calibration bounds required")
    return [(v-ev)*s/(i*(1+ei)), (v+ev)*s/(i*(1-ei))]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nodes", type=int, default=3)
    ap.add_argument("--output", type=Path, default=Path("data/hardware_acquisition.csv"))
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = plan(args.nodes)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    metadata = {"schema": "nonagon.hardware.dc-resistor-plan.v1", "nodes": args.nodes,
        "scalar_configurations": len(rows), "programmable_shunts_required": 2,
        "physical_mapping": "Y_physical = conductance_scale_S * K_normalized",
        "response_mapping": "h_normalized = voltage_V * conductance_scale_S / current_A",
        "loads": "t_normalized = shunt_conductance_S / conductance_scale_S",
        "readout": "one fixed calibrated current source and one fixed voltage channel at node_1",
        "hardware_execution": "not performed; CSV is an acquisition plan",
        "validation_required": ["Known labeled pads; connected passive resistive network; fixed ground",
            "Repeat baseline between intervention groups to bound drift",
            "Measure current, relay path resistance/leakage, shunt conductance and voltmeter calibration",
            "Record settling and repeatability; reserve unused configurations for validation",
            "Propagate all uncertainty; decimal centers are not exact-rational physical certificates"]}
    args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps({"csv": str(args.output), "configurations": len(rows),
                      "hardware_execution": False}))


if __name__ == "__main__":
    main()
