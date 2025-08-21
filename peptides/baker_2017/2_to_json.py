import pandas as pd
from collections import Counter, defaultdict
from Bio import SeqIO
import json
from pathlib import Path
from Bio.PDB import PDBParser, PPBuilder
from tqdm import tqdm

dataset_name = "ha" # botn / ha
dataset_name = "botn" # botn / ha

csv_path = Path(f"{dataset_name}.csv")

df = pd.read_csv(csv_path)
print(f"[Number of unique sequences]: peptide {len(set(df['Sequence']))}, target {len(set(df['target_seq']))}, peptide from pdb {len(set(df['sequence_from_pdb']))}")

seq2categories = defaultdict(list)
peptide_seq2id = dict()
all_target_seqs = set()
for seq, cat, id, target, sequence_from_pdb in zip(df["Sequence"], df["Category"], df["#ID"], df["target_seq"], df["sequence_from_pdb"]):
    if seq not in peptide_seq2id:
        peptide_seq2id[seq] = id
    seq2categories[(seq, target)].append((cat, id, sequence_from_pdb))
    # assert seq == sequence_from_pdb, f"{seq} {sequence_from_pdb}, {sequence_from_pdb in seq}"
    all_target_seqs.add(target)
print(f"[Number of seq pairs]: {len(seq2categories)}")

count = dict(Counter(df["Category"])) # .most_common()
print("[Categories]", count[0], count[1], count[2], count[3], sum(count.values()))

all_target_seqs = list(all_target_seqs)
all_target_seqs.sort()
target_seq2id = dict()
with open(f"{dataset_name}_target_seqs.fasta", "w") as fout:
    for i, seq in enumerate(all_target_seqs):
        fout.write(f">{dataset_name}_target_{i}\n{seq}\n")
        target_seq2id[seq] = f"{dataset_name}_target_{i}"
with open(f"{dataset_name}_binder_seqs.fasta", "w") as fout:
    for seq, id in peptide_seq2id.items():
        fout.write(f">{id}\n{seq}\n")


records = []
conflict_cases = 0
for seq, target in seq2categories:
    peptide_id = peptide_seq2id[seq]
    
    attrs = []
    categories = []
    for cat, id, sequence_from_pdb in seq2categories[(seq, target)]:
        categories.append(cat)
        attrs.append({"category": cat, "id": id, "sequence_from_pdb": sequence_from_pdb})

    categories_set = set(categories)
    if len(categories_set) > 1:
        print(categories_set)
        conflict_cases += 1
        
    records.append({
        "interactors": [
            {'id': target_seq2id[target], 'role': 'target', 'seq': target},
            {'id': peptide_id, 'role': 'binder', 'seq': seq},
        ],
        "attributes": attrs,
        "binding": [int(x) for x in categories][0]
    })

print(f"[number of records]", len(records))
print(f"[number of conflict]", conflict_cases)
    
json.dump(records, open(f"{dataset_name}.json", "w"))
    