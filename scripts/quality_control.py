import os
import hashlib
import pandas as pd
from Bio import SeqIO

MIN_LENGTH = 1400
MAX_N_PCT = 1.0
VALID_BASES = set("ACGTN")

def normalise_sequence(seq_str):
return "".join(seq_str.upper().split())

def sequence_hash(seq_str):
return hashlib.sha256(seq_str.encode()).hexdigest()

def main():
os.makedirs("results", exist_ok=True)
os.makedirs("data/processed", exist_ok=True)

if not os.path.exists("data/raw/candidate_sequences.fasta"):
print("Error: data/raw/candidate_sequences.fasta not found. Run acquire_sequences.py first.")
return

raw_records = list(SeqIO.parse("data/raw/candidate_sequences.fasta", "fasta"))
qc_log = []
qualified_by_species = {}
seen_hashes = set()

print("--- Running Quality Control & Selection ---")
for rec in raw_records:
acc = rec.id
seq = normalise_sequence(str(rec.seq))
length = len(seq)

# Extract species info from description or default
# For simplicity, we parse basic attributes
species = "Unknown"
for word in rec.description.split():
if "OS=" in word:
species = word.replace("OS=", "")

# 1. Length Check
if length < MIN_LENGTH:
qc_log.append({"accession": acc, "status": "EXCLUDE", "reason": "Length < 1400 bp"})
continue

# 2. Invalid Characters Check
invalid = set(seq) - VALID_BASES
if invalid:
qc_log.append({"accession": acc, "status": "EXCLUDE", "reason": f"Invalid characters: {invalid}"})
continue

# 3. Ambiguity Check
n_count = seq.count("N")
n_pct = (n_count / length) * 100
if n_pct > MAX_N_PCT:
qc_log.append({"accession": acc, "status": "EXCLUDE", "reason": f"Excessive N ({n_pct:.2f}%)"})
continue

# 4. Duplicate Check
h = sequence_hash(seq)
if h in seen_hashes:
qc_log.append({"accession": acc, "status": "EXCLUDE", "reason": "Exact sequence duplicate"})
continue
seen_hashes.add(h)

# Passes QC
qc_log.append({"accession": acc, "status": "INCLUDE", "reason": "Passed QC"})

# Group for species selection (placeholder grouping by accession prefix or description)
# In practice, map back to species cleanly

qc_df = pd.DataFrame(qc_log)
qc_df.to_csv("results/dataset_qc_report.csv", index=False)
print("QC report saved to results/dataset_qc_report.csv")

if __name__ == "__main__":
main()


