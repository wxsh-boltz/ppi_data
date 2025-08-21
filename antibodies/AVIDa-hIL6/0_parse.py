import pandas as pd
from collections import defaultdict
import json

path = "AVIDa-hIL6.csv"
dataset_name = "AVIDa-hIL6"
antigen_sequence_path = "antigen_sequences.csv"
saving_path = f"{dataset_name}.json"

# VHH_sequence,label,Ag_label,subject_species,subject_name,subject_sex
df = pd.read_csv(path)
print(len(df), len(set(zip(df["VHH_sequence"], df["Ag_label"]))))
print(sum(df["label"]))

vhh_seq2ids = defaultdict(list)
for i, seq in enumerate(df["VHH_sequence"]):
    vhh_seq2ids[seq].append(i)
print(f"Number of unique VHHs: {len(vhh_seq2ids)}")

vhh_seq2one_id = dict()
with open(f"{dataset_name}_vhh_seqs.fasta", "w") as fout:
    for seq in vhh_seq2ids:
        ids = vhh_seq2ids[seq]
        ids = [f"AVIDa-hIL6_VHH_{idd}" for idd in ids]
        # if len(ids) >= 2:
        #     alias = "|".join(ids[1:])
        # else:
        #     alias = ""
        # fout.write(f">{ids[0]} {alias}\n{seq}\n")
        fout.write(f">{ids[0]}\n{seq}\n")
        vhh_seq2one_id[seq] = ids[0]

df_antigen = pd.read_csv(antigen_sequence_path)
# Ag_label,Ag_sequence
antigen_id2seq = dict(zip(df_antigen["Ag_label"], df_antigen["Ag_sequence"]))
with open(f"{dataset_name}_antigen_seqs.fasta", "w") as fout:
    for antigen_name, antigen_seq in zip(df_antigen["Ag_label"], df_antigen["Ag_sequence"]):
        fout.write(f">{dataset_name}_{antigen_name}\n{antigen_seq}\n")
print(f"Number of antigen: {len(antigen_id2seq)}")
assert len(antigen_id2seq) == len(set(antigen_id2seq.values()))

pair2records = defaultdict(list)
for row in df.to_dict(orient='records'):
    # print(row)
    vhh_seq = row.pop("VHH_sequence")
    pair2records[(vhh_seq, row["Ag_label"])].append(row)


records = []
for vhh_seq, antigen in pair2records:
    interactors = []
    interactors.append({"seq": vhh_seq, "id": vhh_seq2one_id[vhh_seq], "role": "vh"})
    interactors.append({"seq": antigen_id2seq[antigen], "id": f"{dataset_name}_{antigen}", "role": "antigen"})
    binding = set([x["label"] for x in pair2records[vhh_seq, antigen]])
    assert len(binding) == 1
    binding = list(binding)[0]
    if binding == 1:
        binding = True
    else:
        binding = False
    records.append(
        {"interactors": interactors, 
        "attributes": pair2records[vhh_seq, antigen], 
        "binding": binding, 
        "assay_id": dataset_name})

# print(f"Number of records (final)")
print(f"Number of records: {len(records)}")

print("Example:")
print(records[0])

json.dump(records, open(saving_path, "w"))