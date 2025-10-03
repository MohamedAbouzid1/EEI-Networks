#!/usr/bin/env python3
import os

INPUT_FILE = "/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/mm_RNA_seq_data/BRCA_combined_EEIs.txt"
OUTPUT_FILE = "/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main/mm_RNA_seq_data/BRCA_combined_EEIs_unique.txt"


def main() -> None:
    if not os.path.exists(INPUT_FILE):
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    with open(INPUT_FILE, "r") as fin:
        lines = fin.readlines()

    if not lines:
        # Empty input; write empty output
        with open(OUTPUT_FILE, "w") as fout:
            pass
        print(f"No data found. Wrote empty file to: {OUTPUT_FILE}")
        return

    header = lines[0].rstrip("\n")
    col_names = header.split("\t")

    # Assume first two columns are V1 and V2 by name (safer than position-only)
    try:
        v1_idx = col_names.index("V1")
        v2_idx = col_names.index("V2")
    except ValueError:
        # Fallback to first two columns if not labeled
        v1_idx, v2_idx = 0, 1

    seen_pairs = set()
    unique_rows = []

    for raw in lines[1:]:
        row = raw.rstrip("\n")
        if not row:
            continue
        parts = row.split("\t")
        if len(parts) <= max(v1_idx, v2_idx):
            continue
        v1 = parts[v1_idx].strip()
        v2 = parts[v2_idx].strip()
        # Treat interactions as unordered pairs
        pair_key = tuple(sorted((v1, v2)))
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)
        unique_rows.append(row)

    with open(OUTPUT_FILE, "w") as fout:
        fout.write(header + "\n")
        for row in unique_rows:
            fout.write(row + "\n")

    print(f"Input rows (excluding header): {max(0, len(lines)-1):,}")
    print(f"Unique interactions written: {len(unique_rows):,}")
    print(f"Output file: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
