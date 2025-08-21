import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import json
from pathlib import Path
from copy import deepcopy

# [('SARS-CoV2_WT', 12169),
#  ('SARS-CoV1', 3526),
#  ('SARS-CoV2_Omicron-BA1', 2582),
#  ('SARS-CoV2_Omicron-BA2', 2438),
#  ('OC43', 1563),
#  ('HKU1', 1551),
#  ('SARS-CoV2_Omicron-BA2.75', 1459),
#  ('SARS-CoV2_Delta', 1310), # Delta/B.1.617.2
#  ('SARS-CoV2_Beta', 1201),
#  ('SARS-CoV2_Omicron-BA5', 1184),
#  ('SARS-CoV2_Omicron-BA2.12.1', 837),
#  ('Pangolin-GD', 822),
#  ('SARS-CoV2_Omicron-BA1.1', 819),
#  ('SARS-CoV2_Omicron-BA3', 777),
#  ('SARS-CoV2_Alpha', 739),
#  ('SARS-CoV2_Gamma', 738),
#  ('SARS-CoV2_Omicron-BA2.13', 720),
#  ('SARS-CoV2_Omicron-BA4/5', 640),
#  ('SARS-CoV2_Omicron-XBB', 518)]

target_seq_path = "/data/rbg/users/wxsh/esm3/data/antibodies/cov_abdab/ref_seqs/S_variants.fasta"
cov1_target_seq_path = "/data/rbg/users/wxsh/esm3/data/antibodies/cov_abdab/ref_seqs/SARS-CoV1_WT_S.fasta"
raw_data = Path("/data/rbg/users/wxsh/esm3/data/antibodies/cov_abdab/CoV-AbDab_080224.csv")
saving_dir = Path("processed")

# variant_names = ["SARS-CoV2_WT", "SARS-CoV2_Omicron-BA1", "SARS-CoV2_Omicron-BA2"]
variant_names = ["SARS-CoV2_WT", "SARS-CoV2_Omicron-BA1", "SARS-CoV2_Omicron-BA2", "SARS-CoV2_Omicron-BA2.75",\
    "SARS-CoV2_Delta", "SARS-CoV2_Beta", "SARS-CoV2_Omicron-BA5", "SARS-CoV2_Omicron-BA2.12.1", \
    "SARS-CoV2_Alpha", "SARS-CoV2_Gamma", "SARS-CoV1"]

protein_names = ["S; RBD", "S; NTD", "S"] # , ""]

antigen_seqs = dict()
antigen_seq2ids = defaultdict(list)
for record in SeqIO.parse(target_seq_path, "fasta"):
    variant_name = record.id.split("_")[0]
    if len(record.id.split("_")) == 2:
        protein_name = f"S; {record.id.split('_')[1]}"
    else:
        protein_name = "S"
    variant_name = variant_name.replace("Delta/B.1.617.2", "Delta")
    variant_name = variant_name.replace("Beta/B.1.351", "Beta")
    variant_name = variant_name.replace("Gamma/P.1", "Gamma")
    variant_name = variant_name.replace("Alpha/B.1.1.7", "Alpha")

    id = variant_name.replace("/", "-")
    id = f"SARS-CoV2_{id}"
    id = id.replace("BA.1", "BA1")
    id = id.replace("BA.2", "BA2")
    id = id.replace("BA.5", "BA5")
    if id in variant_names:
        antigen_seqs[(id, protein_name)] = str(record.seq)
        # antigen_seq2ids[str(record.seq)].append(record.id)
        antigen_seq2ids[str(record.seq)].append(f"{id}_{protein_name.replace('; ', '_')}")

for record in SeqIO.parse(cov1_target_seq_path, "fasta"):
    variant_name = record.id.split("_")[0]
    if len(record.id.split("_")) == 2:
        protein_name = f"S; {record.id.split('_')[1]}"
    else:
        protein_name = "S"
    id = variant_name.replace("/", "-")
    if id == "WT":
        id = "SARS-CoV1"
    if id in variant_names:
        antigen_seqs[(id, protein_name)] = str(record.seq)
        antigen_seq2ids[str(record.seq)].append(f"{id}_{protein_name.replace('; ', '_')}")

# print(antigen_seqs)
# print(len(antigen_seqs))
# antigen_seq2ids = defaultdict(list)
for key in antigen_seqs:
    print(key, len(antigen_seqs[key]))

#     antigen_seq2ids[antigen_seqs[key]].append(key)
# print(f"Antigen seqL: {len(antigen_seq)}")
# print(len(antigen_seq2ids), len(antigen_seqs))
antigen_seq2id = dict()
for seq in antigen_seq2ids:
    ids = antigen_seq2ids[seq]
    ids.sort()
    antigen_seq2id[seq] = ids[0]
# exit()

def get_variants_list(item):
    items = item.split(";")
    new_items = []
    for x in items:
        new_items.extend([y.strip() for y in x.split(",")])
    return new_items

df = pd.read_csv(raw_data)
df = df.fillna("")

# all_variants = set()
# for _, row in df.iterrows():
#     row = row.to_dict()
#     bind_to = get_variants_list(row["Binds to"])
#     not_bind_to = get_variants_list(row["Doesn't Bind to"])
#     all_variants.update(bind_to)
#     all_variants.update(not_bind_to)
# print(len(all_variants))
# # print(all_variants)
# for v in all_variants:
#     print(v)
# exit()

hchain_seqs = set()
lchain_seqs = set()
records_pos = defaultdict(list)
records_neg = defaultdict(list)
# pair2records = defaultdict(list)
for _, row in df.iterrows():
    row = row.to_dict()
    bind_to = get_variants_list(row["Binds to"])
    bind_to_weak = [strain.replace("(weak)", "").strip() for strain in bind_to if "(weak)" in strain]
    bind_to = [strain for strain in bind_to if "(weak)" not in strain]
    not_bind_to = get_variants_list(row["Doesn't Bind to"])
    
    for variant_name in variant_names:
        if variant_name not in bind_to and variant_name not in not_bind_to and variant_name not in bind_to_weak:
            continue
        if variant_name in bind_to and variant_name in not_bind_to: # don't know which one should I trust
            continue
        if variant_name in bind_to_weak and variant_name in not_bind_to: # don't know which one should I trust
            continue
            
        if variant_name in bind_to:
            binding = 2
        elif variant_name in bind_to_weak:
            binding = 1
        elif variant_name in not_bind_to:
            binding = 0
        else:
            continue

        protein_and_epitope = row["Protein + Epitope"]
        antigen = None
        if protein_and_epitope in protein_names:
            antigen = protein_and_epitope
        else: 
            protein = protein_and_epitope.split(";")[0].strip()
            if protein in protein_names: # using the S protein
                antigen = protein
        if antigen is None:
            continue

        vh = row["VHorVHH"]
        vl = row["VL"]
        if (vh == "ND" or not vh) and (vl == "ND" or not vl): # sometimes there are CDRs?
            continue
        
        attribute = deepcopy(row)
        attribute["binding"] = binding

        if binding:
            records_pos[(vh, vl, antigen_seqs[(variant_name, antigen)])].append(attribute)
        else:
            records_neg[(vh, vl, antigen_seqs[(variant_name, antigen)])].append(attribute)
        if vh != "ND" and vh:
            hchain_seqs.add(vh)
        if vl != "ND" and vl:
            lchain_seqs.add(vl)


hchain_seqs = list(hchain_seqs)
hchain_seqs.sort()
print(max(len(x) for x in hchain_seqs), min(len(x) for x in hchain_seqs))
hchain_seq2id = {seq: f"cov_abdab_vh_{id}" for id, seq in enumerate(hchain_seqs)}
lchain_seqs = list(lchain_seqs)
lchain_seqs.sort()
print(max(len(x) for x in lchain_seqs), min(len(x) for x in lchain_seqs))
lchain_seq2id = {seq: f"cov_abdab_vl_{id}" for id, seq in enumerate(lchain_seqs)}

print(f"antigen_seq2id: {len(antigen_seq2id)}, hchain_seqs: {len(hchain_seqs)}, lchain_seqs: {len(lchain_seqs)}, records_pos: {len(records_pos)}, records_neg: {len(records_neg)}")
# print(len(set((x["VHorVHH"], x["VL"]) for x in rows_filter)))
print(len(set(records_pos.keys()) & set(records_neg.keys()))) # remove this!!!
remove_keys = set(records_pos.keys()) & set(records_neg.keys())
# print(list(remove_keys)[0])
records_pos = {k: v for k, v in records_pos.items() if k not in remove_keys}
records_neg = {k: v for k, v in records_neg.items() if k not in remove_keys}
print(len(set(records_pos.keys()) & set(records_neg.keys()))) # remove this!!!

def add_record(vh, vl, antigen, attributes):
    interactors = []
    if vh != "ND" and vh:
        interactors.append({"id": hchain_seq2id[vh], "seq": vh, "role": "vh"})
    if vl != "ND" and vl:
        interactors.append({"id": lchain_seq2id[vl], "seq": vl, "role": "vl"})
    
    interactors.append({"id": antigen_seq2id[antigen], "seq": antigen, "role": "antigen"})
    assert len(interactors) >= 2, f"{vh}, {vl}, {antigen}"

    attrs_delete_seq = []
    bindings = []
    for attr in attributes:
        attr.pop("VHorVHH", None)
        attr.pop("VL", None)
        attrs_delete_seq.append(attr)
        bindings.append(attr["binding"])

    bindings = list(set(bindings))
    if len(bindings) == 1:
        binding = bindings[0]
    elif 0 in bindings and 1 in bindings or 0 in bindings and 2 in bindings:
        raise ValueError(f"Conflict binding labels: {bindings}")
    elif 1 in bindings and 2 in bindings:
        binding = 2
    else:
        print(bindings)
        exit()

    return {
        "interactors": interactors, 
        "attributes": attrs_delete_seq, 
        "binding": binding,
        "assay_id": "cov_abdab"
        }

outputs = []
for (vh, vl, antigen), attributes in records_pos.items():
    record = add_record(vh, vl, antigen, attributes)
    outputs.append(record)
    # interactors = []
    # if vh != "ND" and vh:
    #     interactors.append({"id": hchain_seq2id[vh], "seq": vh, "role": "vh"})
    # if vl != "ND" and vl:
    #     interactors.append({"id": lchain_seq2id[vl], "seq": vl, "role": "vl"})
    
    # interactors.append({"id": antigen_seq2id[antigen], "seq": antigen, "role": "antigen"})
    # # print(interactors)
    # attrs_delete_seq = []
    # for attr in records_pos[(vh, vl, antigen)]:
    #     attr.pop("VHorVHH", None)
    #     attr.pop("VL", None)
    #     attrs_delete_seq.append(attr)

    # # print(len(attrs_delete_seq))
    # # print(attrs_delete_seq[0])
    # outputs.append({"interactors": interactors, "attributes": attrs_delete_seq, "binding": True})
print(outputs[-1])

for (vh, vl, antigen), attributes in records_neg.items():
    record = add_record(vh, vl, antigen, attributes)
    outputs.append(record)
print(outputs[-1])
# for vh, vl, antigen in records_neg:
#     interactors = []
#     if vh != "ND" and vh:
#         interactors.append({"id": hchain_seq2id[vh], "seq": vh, "role": "vh"})
#     if vl != "ND" and vl:
#         interactors.append({"id": lchain_seq2id[vl], "seq": vl, "role": "vl"})
    
#     interactors.append({"id": antigen_seq2id[antigen], "seq": antigen, "role": "antigen"})
#     attrs_delete_seq = []
#     for attr in records_neg[(vh, vl, antigen)]:
#         attr.pop("VHorVHH", None)
#         attr.pop("VL", None)
#         attrs_delete_seq.append(attr)

#     outputs.append({"interactors": interactors, "attributes": attrs_delete_seq, "binding": False})
print(len(outputs))

json.dump(outputs, open(saving_dir / f"{raw_data.stem}.json", "w"))

with open(saving_dir / f"{raw_data.stem}_antigen_seqs.fasta", "w") as fout:
    for seq in antigen_seq2id:
        fout.write(f">{antigen_seq2id[seq]}\n{seq}\n")

with open(saving_dir / f"{raw_data.stem}_vh_seqs.fasta", "w") as fout:
    for seq in hchain_seq2id:
        fout.write(f">{hchain_seq2id[seq]}\n{seq}\n")
    
with open(saving_dir / f"{raw_data.stem}_vl_seqs.fasta", "w") as fout:
    for seq in lchain_seq2id:
        fout.write(f">{lchain_seq2id[seq]}\n{seq}\n")
    