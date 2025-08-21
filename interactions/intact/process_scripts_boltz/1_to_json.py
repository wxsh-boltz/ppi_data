from Bio import SeqIO
import sys
import pandas as pd
import numpy as np
# import seaborn as sns
from collections import defaultdict, Counter
from Bio import SeqIO
import json
from tqdm import tqdm

input_path = sys.argv[1]

if input_path.endswith(".csv"):
    output_path = input_path.replace(".csv", ".json")
elif input_path.endswith(".txt"):
    output_path = input_path.replace(".txt", ".json")

interaction_attr = [
    "Interaction detection method(s)", 
    "Interaction type(s)", 
    "Source database(s)", 
    # "Interaction identifier(s)",
    "Confidence value(s)",
    "Expansion method(s)",
    "Host organism(s)",
    "Interaction parameter(s)",
    
]
interactor_attr = [
    "Alt. ID(s) interactor",
    "Taxid interactor",
    "Biological role(s) interactor",
    "Experimental role(s) interactor",
    "Type(s) interactor",
    "Stoichiometry(s) interactor",
    # "Feature(s) interactor",
    "Identification method participant"
    
]

def get_uniprotkb_ids(alt_ids):
    return [id.split("uniprotkb:")[1] for id in alt_ids.split("|") if id.startswith("uniprotkb:")]
    for id in alt_ids.split("|"):
        if id.startswith("uniprotkb:"):
            return id.split("uniprotkb:")[1]

# def get_uniprot_seq(uniprot_id, uniprot_id2seq):
#     if uniprot_id in uniprot_id2seq:
#         return uniprot_id2seq[uniprot_id]
#     elif "-PRO_" not in uniprot_id and "-" in uniprot_id:
#         return uniprot_id2seq[uniprot_id.split("-")[0]]
    
def get_seq(ids, id2seq):
    if ids is None:
        return None
    for id in ids:
        if id in id2seq:
            return id2seq[id]

    # for id in ids:
    #     if "-PRO_" not in id and "-" in id: # return the canonical seq
    #         return id2seq[id.split("-")[0]]

def get_intact_ids(alt_ids):
    return [id for id in alt_ids.split("|") if id.startswith("intact:")]
    for id in alt_ids.split("|"):
        if id.startswith("intact:"):
            return id
            # return id.split("intact:")[1]

#### read seqs
uniprot_id2seq = dict()
uniprot_sprot_path = "/data/rbg/users/wxsh/datasets/uniprotkb_20250205/uniprotkb_AND_reviewed_true_2025_04_16.fasta"
# uniprot_sprot_subseq_path = "/data/rbg/users/wxsh/esm3/data/interactions/uniprot/uniprot_sprot_canonical_seqs_and_subseqs.fasta"
print("Reading the sequences from sprot:")
for record in SeqIO.parse(uniprot_sprot_path, "fasta"):
    uniprot_id = record.description.split("|")[1]
    if uniprot_id in uniprot_id2seq:
        assert uniprot_id2seq[uniprot_id] == str(record.seq)
    uniprot_id2seq[uniprot_id] = str(record.seq)
# for record in SeqIO.parse(uniprot_sprot_subseq_path, "fasta"):
#     if "-PRO_" not in record.id:
#         ids = record.id.split("|")
#     else:
#         ids = record.id.split("-PRO_")[0].split("|")
#         pro_id = record.id.split("-PRO_")[1]
#         ids = [x + "-PRO_" + pro_id for x in ids]

#     for idd in ids:
#         uniprot_id2seq[idd] = str(record.seq)

trembl_paths = [
    "/data/rbg/users/wxsh/esm3/data/interactions/intact/uniprot_ids_wo_seq.fasta",
    "/data/rbg/users/wxsh/esm3/data/interactions/intact/uniprot_ids_wo_seq_2.fasta"
]

print("Reading the sequences from trembl (subset):")
for trembl_path in trembl_paths:
    for record in SeqIO.parse(trembl_path, "fasta"):
        if len(record.description.split("|")) >= 2:
            uniprot_id = record.description.split("|")[1]
            # uniprot_id2seq[] = str(record.seq)
        else:
            uniprot_id = record.id
        if uniprot_id in uniprot_id2seq:
            print(f"Mismatch sequence: {uniprot_id}")
            # assert uniprot_id2seq[uniprot_id] == str(record.seq), f"{uniprot_id}\n{str(record.seq)}\n{uniprot_id2seq[uniprot_id]}"
        else:
            uniprot_id2seq[uniprot_id] = str(record.seq)

print("Reading the sequences from IntAct:")
intact_fasta_path = "/data/rbg/users/wxsh/esm3/data/interactions/intact/all_20250613/intact.fasta"
intact_id2seq = dict()
for record in SeqIO.parse(intact_fasta_path, "fasta"):
    intact_id2seq[record.id.replace("INTACT:", "intact:")] = str(record.seq)


#### assign one unique id to each sequence
seq2intact_ids = defaultdict(list)
seq2uniprot_ids = defaultdict(list)
for id, seq in intact_id2seq.items():
    seq2intact_ids[seq].append(id)
for id, seq in uniprot_id2seq.items():
    seq2uniprot_ids[seq].append(id)
seq2unique_id = dict()
for seq in seq2intact_ids:
    ids = seq2intact_ids[seq]
    ids.sort()
    seq2unique_id[seq] = ids[0]
for seq in seq2uniprot_ids:
    if seq in seq2unique_id: # prefer to use the IntAct id
        continue
    ids = seq2uniprot_ids[seq]
    ids.sort()
    seq2unique_id[seq] = ids[0]

if input_path.endswith(".csv"):
    unique_sequence_path = input_path.replace(".csv", "_sequences_unique.fasta")
elif input_path.endswith(".txt"):
    unique_sequence_path = input_path.replace(".txt", "_sequences_unique.fasta")

with open(unique_sequence_path, "w") as fout:
    for seq, id in seq2unique_id.items():
        all_ids = "|".join(seq2intact_ids[seq])
        all_uniprot_ids = "|".join([f"uniprotkb:{id}" for id in seq2uniprot_ids[seq]])
        if all_uniprot_ids:
            all_ids = all_ids + "|" + all_uniprot_ids
        fout.write(f">{id} {all_ids}\n{seq}\n")

# unique_sequence_path = "/data/rbg/users/wxsh/esm3/data/interactions/intact/sequences_unique.fasta"
# seq2unique_id = dict()
# for record in SeqIO.parse(unique_sequence_path, "fasta"):
    # seq2unique_id[str(record.seq)] = record.id

####

if input_path.endswith(".csv"):
    df = pd.read_csv(input_path)
elif input_path.endswith(".txt"):
    df = pd.read_csv(input_path, sep="\t")
print(f"Read {len(df)} lines from {input_path}")

records = defaultdict(list)
fail_to_find_uniprot_id = 0
fail_to_find_seqs = 0
fail_to_find_uniprot_id_set = set()
fail_to_find_seqs_set = set()

for index, row in tqdm(df.iterrows()):
    row_dict = row.to_dict()

    ### get the intact id:
    intact_ids_a = get_intact_ids(row_dict["Alt. ID(s) interactor A"])
    intact_ids_b = get_intact_ids(row_dict["Alt. ID(s) interactor B"])
    
    ### get the uniprot id (if exists):
    uniprotkdb_ids_a = get_uniprotkb_ids(row_dict["#ID(s) interactor A"])
    uniprotkdb_ids_b = get_uniprotkb_ids(row_dict["ID(s) interactor B"])
    uniprotkdb_ids_a.extend(get_uniprotkb_ids(row_dict["Alt. ID(s) interactor A"]))
    uniprotkdb_ids_b.extend(get_uniprotkb_ids(row_dict["Alt. ID(s) interactor B"]))
    
    # if len(uniprotkdb_ids_a) == 0:
    #     uniprotkdb_id_a = 
    # if uniprotkdb_id_b is None:
    #     uniprotkdb_id_b = get_uniprotkb_ids()
    
    ### get seqs:
    # perfer to use the intact seq:
    seq_a = get_seq(intact_ids_a, intact_id2seq)
    seq_b = get_seq(intact_ids_b, intact_id2seq)

    if seq_a is None: # and uniprotkdb_id_a is not None:
        seq_a = get_seq(uniprotkdb_ids_a, uniprot_id2seq)
    if seq_b is None: # and uniprotkdb_id_b is not None:
        seq_b = get_seq(uniprotkdb_ids_b, uniprot_id2seq)
    
    # if uniprotkdb_id_a is None or uniprotkdb_id_b is None:
    #     fail_to_find_uniprot_id += 1
    #     if uniprotkdb_id_a is None:
    #         fail_to_find_uniprot_id_set.add(row_dict["Alt. ID(s) interactor A"])
    #     if uniprotkdb_id_b is None:
    #         fail_to_find_uniprot_id_set.add(row_dict["Alt. ID(s) interactor B"])
    #     continue

    if seq_a is None or seq_b is None:
        if seq_a is None:
            # fail_to_find_seqs_set.update([(x, row_dict["Interaction identifier(s)"]) for x in uniprotkdb_ids_a])
            fail_to_find_seqs_set.update(uniprotkdb_ids_a)
            # fail_to_find_seqs_set.update([(x, row_dict["Interaction identifier(s)"]) for x in intact_ids_a])
        if seq_b is None:
            # fail_to_find_seqs_set.update([(x, row_dict["Interaction identifier(s)"]) for x in uniprotkdb_ids_b])
            fail_to_find_seqs_set.update(uniprotkdb_ids_b)
            # fail_to_find_seqs_set.update([(x, row_dict["Interaction identifier(s)"]) for x in intact_ids_b])
            # fail_to_find_seqs_set.update(intact_ids_b)

        fail_to_find_seqs += 1
        continue
        
    interaction_feats = {}
    for key in interaction_attr:
        interaction_feats[key] = row_dict[key]
        
    interactor_feats_a = {}
    # interactor_feats_a["ids"] = row_dict["Alt. ID(s) interactor A"]
    for key in interactor_attr:
        interactor_feats_a[key] = row_dict[key + " A"]
        
    interactor_feats_b = {}
    # interactor_feats_b["ids"] = row_dict["Alt. ID(s) interactor B"]
    for key in interactor_attr:
        interactor_feats_b[key] = row_dict[key + " B"]

    if seq_a >= seq_b:
        seq_a, seq_b = seq_b, seq_a
        interactor_feats_a, interactor_feats_b = interactor_feats_b, interactor_feats_a
    
    records[(seq_a, seq_b)].append(
        {
            "interaction_attr": interaction_feats, 
            "interactor_attr": [interactor_feats_a, interactor_feats_b]
        }
    )

print(f"fail_to_find_uniprot_id: {fail_to_find_uniprot_id}")
print(f"fail_to_find_seqs: {fail_to_find_seqs}")
with open("fail_to_find_uniprot_id.log", "w") as fout:
    for id in fail_to_find_uniprot_id_set:
        fout.write(f"{id}\n")

with open("fail_to_find_seqs.log", "w") as fout:
    for id in fail_to_find_seqs_set:
        fout.write(f"{id}\n")


records_json = []
for seq_a, seq_b in records:
    id_a = seq2unique_id[seq_a]
    id_b = seq2unique_id[seq_b]
    r = {}
    r["attributes"] = records[(seq_a, seq_b)]
    r["interactors"] = [
        {"id": id_a, "seq": seq_a},
        {"id": id_b, "seq": seq_b}
    ]
    records_json.append(r)

print(f"Number of unique seq pairs {len(records_json)}")

with open(output_path, "w") as fout:
    json.dump(records_json, fout)
