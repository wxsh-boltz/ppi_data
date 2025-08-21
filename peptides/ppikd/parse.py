import csv, json
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from pathlib import Path

dataset_name="ppikd"
path = "PPIKD.csv"
saving_dir = Path("processed")
saving_dir.mkdir(exist_ok=True, parents=True)
df = pd.read_csv(path)

seq_pairs = defaultdict(list)
protein_seqs = set()
peptide_seqs = set()

for index, row in df.iterrows():
    row_dict = row.to_dict()
    protein_seq = row_dict.pop("Protein_Sequence")
    peptide_seq = row_dict.pop("Peptide_Sequence")
    protein_seqs.add(protein_seq)
    peptide_seqs.add(peptide_seq)
    seq_pairs[(protein_seq, peptide_seq)].append(row_dict)
    assert (peptide_seq, protein_seq) not in seq_pairs

protein_seqs = list(protein_seqs)
protein_seqs.sort()
print("protein_seqs", len(protein_seqs), protein_seqs[0])
protein_seq2id = dict()
with open(saving_dir / "target_seqs.fasta", "w") as fout:
    for i, seq in enumerate(protein_seqs):
        id = f"{dataset_name}_target_{i}"
        protein_seq2id[seq] = id
        fout.write(f">{id}\n{seq}\n")

peptide_seqs = list(peptide_seqs)
peptide_seqs.sort()
print("peptide_seqs", len(peptide_seqs), peptide_seqs[0])
peptide_seq2id = dict()
with open(saving_dir / "binder_seqs.fasta", "w") as fout:
    for i, seq in enumerate(peptide_seqs):
        id = f"{dataset_name}_binder_{i}"
        peptide_seq2id[seq] = id
        fout.write(f">{id}\n{seq}\n")

print(f"{len(seq_pairs)} pairs")
records = []
for (protein_seq, peptide_seq), metadata in seq_pairs.items():
    interactors = []
    interactors.append({"id": protein_seq2id[protein_seq], "seq": protein_seq, "role": "target"})
    interactors.append({"id": peptide_seq2id[peptide_seq], "seq": peptide_seq, "role": "binder"})
    records.append({"interactors": interactors, "binding": True, "assay_id": dataset_name, "attributes": metadata})

print("records", len(records))
json.dump(records, open(saving_dir / "data.json", "w"))