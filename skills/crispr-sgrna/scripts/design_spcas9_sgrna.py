#!/usr/bin/env python3
from pathlib import Path
import argparse
import re
import pandas as pd

def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name:
        seqs[name] = "".join(chunks).upper()
    return seqs

def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()

def gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 2)

def homopolymer_warning(seq, n=4):
    for base in "ACGT":
        if base * n in seq:
            return f"homopolymer_{base}{n}+"
    return ""

def score_sgrna(spacer, start, seq_len):
    score = 100
    gc = gc_content(spacer)
    warnings = []

    if gc < 40:
        score -= 20
        warnings.append("low_GC")
    elif gc > 70:
        score -= 20
        warnings.append("high_GC")

    hp = homopolymer_warning(spacer)
    if hp:
        score -= 10
        warnings.append(hp)

    if "N" in spacer:
        score -= 50
        warnings.append("ambiguous_base")

    # Prefer targets in the first half for knockout from CDS sequence
    relative_pos = start / max(seq_len, 1)
    if relative_pos <= 0.35:
        score += 10
    elif relative_pos >= 0.75:
        score -= 10
        warnings.append("late_position")

    # Prefer G near 5' for transcription compatibility, weak bonus
    if spacer.startswith("G"):
        score += 3

    return score, ";".join(warnings) if warnings else "OK"

def find_spcas9_targets(seq_name, seq):
    records = []
    seq_len = len(seq)

    # plus strand: N20 NGG
    for m in re.finditer(r"(?=([ACGTN]{20}([ACGT]GG)))", seq):
        spacer = m.group(1)[:20]
        pam = m.group(2)
        start = m.start(1) + 1
        end = start + 22
        score, warn = score_sgrna(spacer, start, seq_len)
        records.append({
            "seq_id": seq_name,
            "strand": "+",
            "start_1based": start,
            "end_1based": end,
            "spacer_20nt": spacer,
            "pam": pam,
            "gc_percent": gc_content(spacer),
            "score": score,
            "warning": warn
        })

    # minus strand: CCN on plus sequence corresponds to NGG PAM on reverse strand
    for m in re.finditer(r"(?=((CC[ACGT])[ACGTN]{20}))", seq):
        protospacer_plus = m.group(1)[3:23]
        spacer = revcomp(protospacer_plus)
        pam = revcomp(m.group(2))
        start = m.start(1) + 1
        end = start + 22
        score, warn = score_sgrna(spacer, start, seq_len)
        records.append({
            "seq_id": seq_name,
            "strand": "-",
            "start_1based": start,
            "end_1based": end,
            "spacer_20nt": spacer,
            "pam": pam,
            "gc_percent": gc_content(spacer),
            "score": score,
            "warning": warn
        })

    return records

def main():
    parser = argparse.ArgumentParser(description="Design SpCas9 NGG sgRNAs from FASTA.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output CSV")
    parser.add_argument("--top", type=int, default=10, help="Number of top targets to print")
    args = parser.parse_args()

    seqs = read_fasta(args.input)
    all_records = []
    for name, seq in seqs.items():
        all_records.extend(find_spcas9_targets(name, seq))

    df = pd.DataFrame(all_records)
    if df.empty:
        print("No SpCas9 NGG sgRNA targets found.")
        df.to_csv(args.output, index=False)
        return

    df = df.sort_values(["score", "gc_percent"], ascending=[False, True])
    df["rank"] = range(1, len(df) + 1)
    cols = ["rank", "seq_id", "strand", "start_1based", "end_1based",
            "spacer_20nt", "pam", "gc_percent", "score", "warning"]
    df = df[cols]
    df.to_csv(args.output, index=False)

    print(f"Found {len(df)} candidate sgRNAs.")
    print(df.head(args.top).to_string(index=False))

if __name__ == "__main__":
    main()
