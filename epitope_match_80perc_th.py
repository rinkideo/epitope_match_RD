__version__ = "1.0.0"
import pandas as pd
import gzip
import argparse
from Bio.Seq import Seq
from collections import defaultdict
from time import time
import os
import warnings
from Bio import BiopythonWarning
from openpyxl import Workbook

warnings.simplefilter('ignore', BiopythonWarning)

# === Reverse complement helper ===
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

# === FASTQ reader ===
def read_fastq_sequences(file_path):
    with gzip.open(file_path, "rt") as handle:
        while True:
            try:
                _ = next(handle)
                seq = next(handle).strip()
                _ = next(handle)
                _ = next(handle)
                yield seq
            except StopIteration:
                break

# === Load HXB2 reference ===
def load_reference(fasta_file):
    with open(fasta_file) as f:
        lines = f.readlines()
    return ''.join(line.strip() for line in lines if not line.startswith(">"))\
.upper()

# === Translate and compare AA sequences ===
def compare_aa_variants(ref_nt_seq, matched_nt_windows):
    ref_aa = str(Seq(ref_nt_seq[:len(ref_nt_seq)//3*3]).translate())
    variant_counter = defaultdict(list)

    for nt_window, identity in matched_nt_windows:
        nt_trimmed = nt_window[:len(nt_window)//3*3]
        aa_window = str(Seq(nt_trimmed).translate())
        variant = ''.join('.' if a == b else b for a, b in zip(ref_aa, aa_window))
        variant_counter[variant].append(identity)

    variant_rows = []
    for variant, identities in sorted(variant_counter.items(), key=lambda x: -len(x[1])):
        count = len(identities)
        percent = round(100 * count / sum(len(v) for v in variant_counter.values()), 2)
        min_identity = round(min(identities) * 100, 2)
        avg_identity = round(sum(identities) / len(identities) * 100, 2)
        row = list(variant)
        variant_rows.append((row, percent, count, min_identity, avg_identity))

    return list(ref_aa), variant_rows

# === Process a single epitope ===
def process_epitope(ref_seq, reads, epitope_name):
    matched_windows = []
    epitope_len = len(ref_seq)
    for read in reads:
        if len(read) >= epitope_len:
            for i in range(len(read) - epitope_len + 1):
                window = read[i:i + epitope_len]
                identity = sum(a == b for a, b in zip(window, ref_seq)) / epitope_len
                if identity >= 0.65:
                    matched_windows.append((window, identity))

    if not matched_windows:
        return None, 0

    ref_aa, variant_rows = compare_aa_variants(ref_seq, matched_windows)
    table = []
    for i, (variant, percent, count, min_identity, avg_identity) in enumerate(variant_rows):
        sample = epitope_name if i == 0 else ""
        table.append([sample] + variant + [percent, count, min_identity, avg_identity])
    header = ["Epitope"] + ref_aa + ["Percent", "Reads", "MinIdentity", "AvgIdentity"]
    return pd.DataFrame(table, columns=header), len(matched_windows)

# === Main logic ===
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--epitopes", required=True)
    parser.add_argument("--output-dir", default="./output")
    args = parser.parse_args()

    sample_id = args.sample_id
    r1_path = args.r1
    r2_path = args.r2
    fasta_file = args.ref
    epitope_excel = args.epitopes
    outdir = args.output_dir

    os.makedirs(outdir, exist_ok=True)

    hxb2 = load_reference(fasta_file)
    r1_reads = list(read_fastq_sequences(r1_path))
    r2_reads = list(read_fastq_sequences(r2_path))
    all_reads = r1_reads + [reverse_complement(r) for r in r2_reads]

    epitopes = pd.read_excel(epitope_excel)
    epitopes = epitopes.dropna(subset=["CLADEB_GENOME_START", "CLADEB_GENOME_END"])
    epitopes["CLADEB_GENOME_START"] = epitopes["CLADEB_GENOME_START"].astype(int)
    epitopes["CLADEB_GENOME_END"] = epitopes["CLADEB_GENOME_END"].astype(int)

    output_tables = []
    unmatched = []
    summary = []

    print(f"⏳ Starting processing for sample: {sample_id}")
    for _, row in epitopes.iterrows():
        name = row["Fernando_Approved_Name"]
        start, end = row["CLADEB_GENOME_START"], row["CLADEB_GENOME_END"]
        ref_nt_seq = hxb2[start - 1:end]

        epitope_start = time()
        df_epitope, count = process_epitope(ref_nt_seq, all_reads, name)
        print(f"🧬 {sample_id} | {name}: {count} matches found in {round(time() - epitope_start, 2)} sec")

        if df_epitope is not None:
            df_epitope.insert(0, "SampleID", sample_id)
            output_tables.append((name, df_epitope))
            summary.append({"SampleID": sample_id, "Epitope": name, "MatchedReads": count})
        else:
            unmatched.append({"SampleID": sample_id, "Epitope": name})
            summary.append({"SampleID": sample_id, "Epitope": name, "MatchedReads": 0})

    # Write outputs
    if output_tables:
        out_path = os.path.join(outdir, f"{sample_id}_variants_{sample_id}.xlsx")
        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            for i, (epitope, df) in enumerate(output_tables):
                sheetname = (epitope[:25] + f"_{i}")[:31]
                df.to_excel(writer, index=False, sheet_name=sheetname)
        print(f"✅ Variants written to: {out_path}")

# === Run main ===
if __name__ == "__main__":
    main()
