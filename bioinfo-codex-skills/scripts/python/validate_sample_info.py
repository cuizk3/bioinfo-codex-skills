#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd

def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: validate_sample_info.py sample_info.csv")

    path = Path(sys.argv[1])
    df = pd.read_csv(path)

    required = {"sample", "group", "fastq1", "fastq2"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    if df["sample"].duplicated().any():
        dup = df.loc[df["sample"].duplicated(), "sample"].tolist()
        raise SystemExit(f"Duplicated sample names: {dup}")

    group_counts = df.groupby("group")["sample"].count()
    small = group_counts[group_counts < 2]
    if len(small) > 0:
        print(f"Warning: groups with fewer than 2 replicates: {small.to_dict()}", file=sys.stderr)

    missing_files = []
    for col in ["fastq1", "fastq2"]:
        for f in df[col]:
            if not Path(f).exists():
                missing_files.append(f)
    if missing_files:
        print("Warning: these FASTQ files do not exist yet:", file=sys.stderr)
        for f in missing_files[:20]:
            print(f"  {f}", file=sys.stderr)
        if len(missing_files) > 20:
            print(f"  ... and {len(missing_files) - 20} more", file=sys.stderr)

    print("sample_info validation finished.")

if __name__ == "__main__":
    main()
