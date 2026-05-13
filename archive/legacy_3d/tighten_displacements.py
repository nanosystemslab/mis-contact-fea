#!/usr/bin/env python3
"""Compute tightened INITIAL_GAP and DISP per shape based on the
previous run's force_displacement.csv (contact onset and release).

Rules:
  cap     : onset target = --cap-onset    (default 90)
  capcone : onset target = --capcone-onset (default 130)
  others  : new onset = current_onset - 5  (start 5 µm earlier)
            new disp  = new_onset + window + 5  (stop 5 µm after release)

For shapes given an explicit onset target, the new INITIAL_GAP is shifted
so that contact engages at the target disp, and DISP is sized to cover
window + 5 µm post-release margin.

The relationship:
   onset_disp = initial_gap + (z_apex_bot - z_contact)
where (z_apex_bot - z_contact) is a shape constant. So
   new_initial_gap = current_initial_gap + (new_onset_disp - current_onset_disp)

Outputs env-var lines that can be eval'd:
   export INITIAL_GAP_<shape>=<value>
   export DISP_<shape>=<value>
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def load_force_csv(path):
    d, f = [], []
    with path.open() as fh:
        reader = csv.reader(fh)
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                try:
                    d.append(float(row[0])); f.append(float(row[1]))
                except ValueError:
                    pass
    return d, f


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("runs_all_shapes"))
    parser.add_argument(
        "--shapes", nargs="+",
        default=["sphere", "oblate", "prolate", "cone", "cap", "capcone", "doublecone"],
    )
    parser.add_argument("--contact-threshold", type=float, default=1.0)
    parser.add_argument("--pre-margin", type=float, default=2.0,
                        help="Absolute pre-contact run-up in µm (= new onset disp). "
                             "Smaller = less wasted free-travel before contact. "
                             "Use ~1–5 µm to give Newton a brief warm-up before contact.")
    parser.add_argument("--post-margin", type=float, default=5.0,
                        help="µm after release in the new DISP.")
    parser.add_argument("--cap-onset", type=float, default=None,
                        help="Override new onset disp specifically for cap "
                             "(default: same --pre-margin as everyone else).")
    parser.add_argument("--capcone-onset", type=float, default=None,
                        help="Override new onset disp specifically for capcone "
                             "(default: same --pre-margin as everyone else).")
    parser.add_argument("--export", choices=["env", "table"], default="table")
    args = parser.parse_args()

    targets = {}
    if args.cap_onset is not None:
        targets["cap"] = args.cap_onset
    if args.capcone_onset is not None:
        targets["capcone"] = args.capcone_onset

    if args.export == "table":
        print(f"{'shape':<12} {'cur_gap':>10} {'cur_disp':>10} {'cur_onset':>11} "
              f"{'cur_release':>13} {'window':>9}  → "
              f"{'new_onset':>10} {'new_gap':>10} {'new_disp':>10}")
        print("-" * 110)

    env_lines = []
    for shape in args.shapes:
        meta_path = args.root / shape / "mesh" / "mesh_tags.json"
        csv_path = args.root / shape / "force_displacement.csv"
        if not meta_path.exists() or not csv_path.exists():
            if args.export == "table":
                print(f"{shape:<12}  (missing data)")
            continue

        meta = json.loads(meta_path.read_text())
        current_gap = float(meta.get("initial_gap_um", 0.0))

        d, f = load_force_csv(csv_path)
        if not d:
            continue
        current_disp = max(d)

        onset_idx = next(
            (i for i, v in enumerate(f) if abs(v) >= args.contact_threshold), None,
        )
        release_idx = None
        if onset_idx is not None:
            for i in range(len(f) - 1, -1, -1):
                if abs(f[i]) >= args.contact_threshold:
                    release_idx = i
                    break
        if onset_idx is None or release_idx is None:
            if args.export == "table":
                print(f"{shape:<12}  (no contact detected)")
            continue

        current_onset = d[onset_idx]
        current_release = d[release_idx]
        window = current_release - current_onset

        # New onset = absolute pre-contact run-up (default 2 µm), or
        # an explicit per-shape override.
        if shape in targets:
            new_onset = targets[shape]
        else:
            new_onset = args.pre_margin
        new_gap = current_gap + (new_onset - current_onset)
        new_disp = (new_onset - new_gap + current_gap)  # placeholder, redone next line
        new_disp = (new_onset + window + args.post_margin) - new_gap + current_gap
        # Simpler: total descent needed = new_onset - new_gap + current_gap + window + post_margin
        # i.e., new_disp covers from start of new sim (which is at new_gap) past onset/release.
        # Easier: new_disp = (new_onset + window + post_margin) - (new_gap - current_gap)
        # since the initial-gap shift moves the "0" of the sim.
        # The sim still starts at disp=0 and runs until disp=new_disp. We need
        # new_disp such that current_release shifted by (new_gap - current_gap)
        # plus post_margin is reached.
        # current_release in the NEW frame = current_release - (current_gap - new_gap)
        # ie current_release + (new_gap - current_gap) reads the same dist along
        # the new sim's axis... actually new_onset = current_onset + (new_gap - current_gap)
        # so new_release = current_release + (new_gap - current_gap) = new_onset + window.
        # We want sim to end at new_release + post_margin.
        new_disp = new_onset + window + args.post_margin

        env_lines.append((shape, new_gap, new_disp))

        if args.export == "table":
            print(f"{shape:<12} {current_gap:>10.2f} {current_disp:>10.2f} "
                  f"{current_onset:>11.2f} {current_release:>13.2f} {window:>9.2f}  → "
                  f"{new_onset:>10.2f} {new_gap:>10.2f} {new_disp:>10.2f}")

    if args.export == "env":
        for shape, gap, disp in env_lines:
            print(f"INITIAL_GAP_{shape}={gap:.3f}")
            print(f"DISP_{shape}={disp:.3f}")


if __name__ == "__main__":
    main()
