import pandas as pd
import json
from collections import defaultdict
from Bio import SeqIO

def get_seq2id(path):
    seq2id = dict()
    for record in SeqIO.parse(path, "fasta"):
        seq2id[str(record.seq)] = record.id
    return seq2id

raw_data_path = "vincent_data_clean.csv"
df = pd.read_csv(raw_data_path)
df = df.fillna("")

hchain_seq2id = get_seq2id("heavy_sequences.fasta")
lchain_seq2id = get_seq2id("light_sequences.fasta")
antigen_seq2id = get_seq2id("antigen_sequences.fasta")


records = defaultdict(list)
for _, row in df.iterrows():
    row = row.to_dict()
    hchain = row.pop("heavy_sequence")
    lchain = row.pop("light_sequence")
    antigen = row.pop("gpt4output_sequence_from_accession")
    if (hchain or lchain) and antigen:
        records[(hchain, lchain, antigen)].append(row)
print(len(df), len(records))

output = []
for (hchain, lchain, antigen), attributes in records.items():
    interactors = []
    if hchain:
        interactors.append({"seq": hchain, "id": hchain_seq2id[hchain], "role": "hchain"})
    if lchain:
        interactors.append({"seq": lchain, "id": lchain_seq2id[lchain], "role": "lchain"})
    if antigen:
        interactors.append({"seq": antigen, "id": antigen_seq2id[antigen], "role": "antigen"})
    output.append({"interactors": interactors, "attributes": attributes})

print(len(output))
json.dump(output, open(raw_data_path.replace(".csv", ".json"), "w"))