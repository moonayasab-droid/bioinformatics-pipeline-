import os
import pandas as pd
from Bio import Entrez, SeqIO

# IMPORTANT: Replace with your actual email address for NCBI Entrez
Entrez.email = "skyhighsab@gmail.com"
Entrez.tool = "comparative_16s_pipeline"

TARGET_TAXA = [
{"phylum": "Pseudomonadota", "family": "Enterobacteriaceae", "genus": "Escherichia", "species": "Escherichia coli"},
{"phylum": "Pseudomonadota", "family": "Enterobacteriaceae", "genus": "Escherichia", "species": "Escherichia fergusonii"},
{"phylum": "Pseudomonadota", "family": "Enterobacteriaceae", "genus": "Klebsiella", "species": "Klebsiella pneumoniae"},
{"phylum": "Pseudomonadota", "family": "Enterobacteriaceae", "genus": "Klebsiella", "species": "Klebsiella oxytoca"},
{"phylum": "Pseudomonadota", "family": "Enterobacteriaceae", "genus": "Klebsiella", "species": "Klebsiella aerogenes"},
{"phylum": "Pseudomonadota", "family": "Enterobacteriaceae", "genus": "Salmonella", "species": "Salmonella enterica"},
{"phylum": "Pseudomonadota", "family": "Enterobacteriaceae", "genus": "Salmonella", "species": "Salmonella bongori"},
{"phylum": "Pseudomonadota", "family": "Pseudomonadaceae", "genus": "Pseudomonas", "species": "Pseudomonas aeruginosa"},
{"phylum": "Pseudomonadota", "family": "Pseudomonadaceae", "genus": "Pseudomonas", "species": "Pseudomonas putida"},
{"phylum": "Pseudomonadota", "family": "Pseudomonadaceae", "genus": "Pseudomonas", "species": "Pseudomonas fluorescens"},
{"phylum": "Pseudomonadota", "family": "Pseudomonadaceae", "genus": "Pseudomonas", "species": "Pseudomonas syringae"},
{"phylum": "Bacillota", "family": "Bacillaceae", "genus": "Bacillus", "species": "Bacillus subtilis"},
{"phylum": "Bacillota", "family": "Bacillaceae", "genus": "Bacillus", "species": "Bacillus cereus"},
{"phylum": "Bacillota", "family": "Bacillaceae", "genus": "Bacillus", "species": "Bacillus anthracis"},
{"phylum": "Bacillota", "family": "Bacillaceae", "genus": "Bacillus", "species": "Bacillus amyloliquefaciens"},
{"phylum": "Bacillota", "family": "Staphylococcaceae", "genus": "Staphylococcus", "species": "Staphylococcus aureus"},
{"phylum": "Bacillota", "family": "Staphylococcaceae", "genus": "Staphylococcus", "species": "Staphylococcus epidermidis"},
{"phylum": "Bacillota", "family": "Staphylococcaceae", "genus": "Staphylococcus", "species": "Staphylococcus lugdunensis"}
]

def main():
os.makedirs("data/raw", exist_ok=True)
all_records = []
metadata_rows = []

print("--- Starting NCBI Candidate Acquisition ---")
for taxon in TARGET_TAXA:
sp = taxon["species"]
# Strict server-side query filters
query = f'"{sp}"[ORGN] AND "16S ribosomal RNA"[TITL] AND biomol_rRNA[PROP] AND 1400:2000[SLEN]'
print(f"Querying: {sp}")

try:
handle = Entrez.esearch(db="nuccore", term=query, retmax=20)
search_results = Entrez.read(handle)
handle.close()

id_list = search_results.get("IdList", [])
if not id_list:
continue

# Fetch raw records
fetch_handle = Entrez.efetch(db="nuccore", id=id_list, rettype="gb", retmode="text")
records = list(SeqIO.parse(fetch_handle, "genbank"))
fetch_handle.close()

for rec in records:
acc = rec.id
seq_str = str(rec.seq)
all_records.append(rec)

metadata_rows.append({
"sequence_id": f"{taxon['genus']}_{taxon['species'].split()[-1]}_{acc}",
"accession": acc,
"organism": sp,
"species": sp,
"genus": taxon["genus"],
"family": taxon["family"],
"order": "NA",
"class": "NA",
"phylum": taxon["phylum"],
"sequence_length": len(seq_str),
"source_database": "NCBI",
"retrieval_date": "2026-03-30"
})
except Exception as e:
print(f"Error processing {sp}: {e}")

# Save raw outputs
SeqIO.write(all_records, "data/raw/candidate_sequences.fasta", "fasta")
meta_df = pd.DataFrame(metadata_rows)
meta_df.to_csv("data/raw/candidate_metadata.csv", index=False)
print(f"Acquisition complete. Raw records saved to data/raw/")

if __name__ == "__main__":
main()


