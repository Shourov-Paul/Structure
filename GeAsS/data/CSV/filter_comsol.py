#!/usr/bin/env python3
"""
COMSOL CSV Filter and Branch Tracker
Filters a COMSOL-generated CSV file by a value range in the 4th column,
then selects exactly one row per index (column 1) by tracking the physical
branch where column 4 is ascending and column 3 is descending.
Supports carriage return (\r) and other line endings automatically.
"""

import os
import csv
import argparse

def detect_line_ending(file_path):
    with open(file_path, 'rb') as f:
        chunk = f.read(4096)
        if b'\r\n' in chunk:
            return '\r\n'
        elif b'\r' in chunk:
            return '\r'
        else:
            return '\n'

def parse_csv(file_path, newline):
    with open(file_path, 'r', newline=newline, encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f)
        return [row for row in reader if row]

def filter_and_track(rows, min_val, max_val):
    """
    Groups rows by column 1 (index).
    - First index: uses [min_val, max_val] to identify the starting branch,
      picks the candidate with the lowest col4 in that range.
    - All subsequent indices: tracks the physical branch freely (no range
      restriction) by preferring candidates where col4 is ascending AND col3
      is descending, choosing the one with the smallest combined normalised
      distance in both col3 and col4 from the previous tracked values.
    """
    grouped = {}
    for row in rows:
        if len(row) < 4:
            continue
        try:
            idx  = float(row[0])
            col3 = float(row[2])
            col4 = float(row[3])
            grouped.setdefault(idx, []).append((col3, col4, row))
        except ValueError:
            continue

    if not grouped:
        return []

    sorted_indices = sorted(grouped.keys())

    tracked = []
    history = []  # (col3, col4) of each tracked point for velocity extrapolation
    tol = 0.01    # 1% tolerance on col4 to handle near-flat mode crossings

    for idx in sorted_indices:
        candidates = grouped[idx]

        if not history:
            # First index: use range to identify starting branch
            in_range = [c for c in candidates if min_val <= c[1] <= max_val]
            if not in_range:
                continue
            best = min(in_range, key=lambda c: c[1])

        else:
            prev_col3, prev_col4 = history[-1]

            # Extrapolate expected next values from recent velocity
            if len(history) >= 2:
                exp_col3 = prev_col3 + (history[-1][0] - history[-2][0])
                exp_col4 = prev_col4 + (history[-1][1] - history[-2][1])
            else:
                exp_col3, exp_col4 = prev_col3, prev_col4

            # Candidates satisfying both physical trends (no range restriction)
            valid = [c for c in candidates
                     if c[1] > prev_col4 * (1 - tol) and c[0] < prev_col3]
            if valid:
                # Pick closest to extrapolated position
                def score(c, ec3=exp_col3, ec4=exp_col4):
                    d3 = abs(c[0] - ec3) / (abs(ec3) + 1e-30)
                    d4 = abs(c[1] - ec4) / (abs(ec4) + 1e-30)
                    return d3 + d4
                best = min(valid, key=score)
            else:
                # Fallback: nearest neighbour in col4
                best = min(candidates, key=lambda c: abs(c[1] - prev_col4))

        tracked.append(best[2])
        history.append((best[0], best[1]))

    return tracked

def main():
    parser = argparse.ArgumentParser(
        description="Filter and track physical branch in COMSOL CSV files.")
    parser.add_argument("-i", "--input",  help="Path to the input CSV file")
    parser.add_argument("-o", "--output", help="Path to the output CSV file (optional)")
    parser.add_argument("--min", type=float, help="Minimum value for the 4th column")
    parser.add_argument("--max", type=float, help="Maximum value for the 4th column")
    args = parser.parse_args()

    if not args.input:
        print("=== COMSOL CSV Filter & Branch Tracker ===")
        args.input = input("Please enter the file path: ").strip().strip('"\'')

    if not os.path.exists(args.input):
        print(f"Error: File '{args.input}' does not exist.")
        return

    if args.min is None:
        args.min = float(input("Enter minimum value for the 4th column (e.g. 6.68e-12): ").strip())

    if args.max is None:
        args.max = float(input("Enter maximum value for the 4th column (e.g. 1.71e-11): ").strip())

    if args.output is None:
        default_out = os.path.splitext(args.input)[0] + "_filtered.csv"
        out_input = input(f"Enter output CSV path [default: {default_out}]: ").strip()
        args.output = out_input if out_input else default_out

    # 1. Detect line endings
    newline = detect_line_ending(args.input)
    print(f"\n[1/3] Detected line ending: {repr(newline)}")

    # 2. Parse CSV
    print("[2/3] Reading input CSV...")
    rows = parse_csv(args.input, newline)
    print(f"      Read {len(rows)} rows.")

    # 3. Filter range + track physical branch (ascending col4, descending col3)
    print(f"[3/3] Filtering col4 [{args.min:.3e}, {args.max:.3e}] and tracking "
          f"physical branch (col4 asc, col3 desc)...")
    tracked = filter_and_track(rows, args.min, args.max)
    print(f"      Tracked {len(tracked)} rows (one per index).")

    # Write output
    print(f"\nWriting results to: '{args.output}'")
    with open(args.output, 'w', newline='') as outfile:
        writer = csv.writer(outfile, lineterminator=newline)
        writer.writerows(tracked)

    print("Done.\n")

if __name__ == "__main__":
    main()
