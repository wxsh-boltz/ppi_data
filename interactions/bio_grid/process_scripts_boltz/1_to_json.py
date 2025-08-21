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

# be consistent with intact
interaction_attr = {
    "Interaction Detection Method": "Interaction detection method(s)", 
    "Interaction Types": "Interaction type(s)", 
    "Source Database": "Source database(s)", 
    "Confidence Values": "Confidence value(s)", 
}
interactor_attr = {
    "Taxid Interactor": "Taxid interactor",
}

def get_uniprot_ids(alt_id):
    alt_id_split = alt_id.split("|")
    uniprot_ids = [x.replace("uniprot/swiss-prot:", "") for x in alt_id_split if "uniprot/swiss-prot:" in x]
    ncbi_ids = [x.replace("refseq:", "") for x in alt_id_split if "refseq:" in x]
    # if len(uniprot_ids) > 0:
    #     uniprot_ids = [uniprot_id.replace("uniprot/swiss-prot:", "") for uniprot_id in uniprot_ids]
    #     # uniprot_id = uniprot_ids[0].replace("uniprot/swiss-prot:", "")
    # else:
    #     uniprot_id = None
    # if len(ncbi_ids) > 0:
    #     ncbi_id = ncbi_ids[0].replace("refseq:", "")
    # else:
    #     ncbi_id = None
    return uniprot_ids, ncbi_ids

def get_seq(uniprot_ids, ncbi_ids, uniprot_id2seq, ncbi_id2seq):
    # if uniprot_id is not None:
    for uniprot_id in uniprot_ids:
        if uniprot_id in uniprot_id2seq:
            return uniprot_id2seq[uniprot_id]
        elif "-PRO_" not in uniprot_id and "-" in uniprot_id: # use the canonical
            return uniprot_id2seq[uniprot_id.split("-")[0]]
    
    for ncbi_id in ncbi_ids:
        if ncbi_id in ncbi_id2seq:
            return ncbi_id2seq[ncbi_id]

    # if ncbi_id is not None and ncbi_id in ncbi_id2seq:
    #     return ncbi_id2seq[ncbi_id]


#### read seqs
uniprot_id2seq = dict()
uniprot_sprot_path = "/data/rbg/users/wxsh/datasets/uniprotkb_20250205/uniprotkb_AND_reviewed_true_2025_04_16.fasta"
# uniprot_sprot_subseq_path = "/data/rbg/users/wxsh/esm3/data/interactions/uniprot/uniprot_sprot_canonical_seqs_and_subseqs.fasta"
print("Reading the sequences from sprot:")
for record in SeqIO.parse(uniprot_sprot_path, "fasta"):
    uniprot_id2seq[record.description.split("|")[1]] = str(record.seq)
# for record in SeqIO.parse(uniprot_sprot_subseq_path, "fasta"):
#     if "-PRO_" not in record.id:
#         ids = record.id.split("|")
#     else:
#         ids = record.id.split("-PRO_")[0].split("|")
#         pro_id = record.id.split("-PRO_")[1]
#         ids = [x + "-PRO_" + pro_id for x in ids]

#     for idd in ids:
#         uniprot_id2seq[idd] = str(record.seq)
uniprot_extended_paths = [
    "/data/rbg/users/wxsh/esm3/data/interactions/bio_grid/process_scripts_boltz/fail_to_find_uniprot_id.fasta",
]
for path in uniprot_extended_paths:
    for record in SeqIO.parse(path, "fasta"):
        uniprot_id = record.id
        if uniprot_id not in uniprot_id2seq:
            uniprot_id2seq[uniprot_id] = str(record.seq)


ncbi_path = "/data/rbg/users/wxsh/esm3/data/interactions/bio_grid/process_scripts_boltz/refseq_ids.fasta"
print("Reading the sequences from ncbi (subset):")
ncbi_id2seq = dict()
for record in SeqIO.parse(ncbi_path, "fasta"):
    if record.id.split(".")[0] in ncbi_id2seq:
        assert ncbi_id2seq[record.id.split(".")[0]] == str(record.seq)
    ncbi_id2seq[record.id.split(".")[0]] = str(record.seq)

# unique_sequence_path = "/data/rbg/users/wxsh/esm3/data/interactions/bio_grid/process_scripts_boltz/sequences_unique.fasta"
# seq2unique_id = dict()
# for record in SeqIO.parse(unique_sequence_path, "fasta"):
#     seq2unique_id[str(record.seq)] = record.id


#### assign one unique id to each sequence
seq2ncbi_ids = defaultdict(list)
seq2uniprot_ids = defaultdict(list)
for id, seq in ncbi_id2seq.items():
    seq2ncbi_ids[seq].append(id)
for id, seq in uniprot_id2seq.items():
    seq2uniprot_ids[seq].append(id)
seq2unique_id = dict()
for seq in seq2uniprot_ids:
    ids = seq2uniprot_ids[seq]
    ids.sort()
    seq2unique_id[seq] = ids[0]
for seq in seq2ncbi_ids:
    if seq in seq2unique_id: # prefer to use the UniProt id
        continue
    ids = seq2ncbi_ids[seq]
    ids.sort()
    seq2unique_id[seq] = ids[0]

if input_path.endswith(".csv"):
    unique_sequence_path = input_path.replace(".csv", "_sequences_unique.fasta")
elif input_path.endswith(".txt"):
    unique_sequence_path = input_path.replace(".txt", "_sequences_unique.fasta")

with open(unique_sequence_path, "w") as fout:
    for seq, id in seq2unique_id.items():
        all_ids = "|".join([f"refseq:{id}" for id in seq2ncbi_ids[seq]])
        all_uniprot_ids = "|".join([f"uniprotkb:{id}" for id in seq2uniprot_ids[seq]])
        if all_uniprot_ids:
            all_ids = all_ids + "|" + all_uniprot_ids
        fout.write(f">{id} {all_ids}\n{seq}\n")

####
if input_path.endswith(".csv"):
    df = pd.read_csv(input_path)
elif input_path.endswith(".txt"):
    df = pd.read_csv(input_path, sep="\t")
print(f"Read {len(df)} lines from {input_path}")

records = defaultdict(list)
fail_to_find_uniprot_id = 0
fail_to_find_seqs_cnt = 0
fail_to_find_seqs_uniprot = set()
fail_to_find_seqs_ncbi = set()

fail_to_find_uniprot_id_set = set()

for index, row in tqdm(df.iterrows()):
    row_dict = row.to_dict()

    uniprot_ids_a, ncbi_ids_a = get_uniprot_ids(row_dict["Alt IDs Interactor A"])
    uniprot_ids_b, ncbi_ids_b = get_uniprot_ids(row_dict["Alt IDs Interactor B"])

    if (len(uniprot_ids_a) == 0 and len(ncbi_ids_a) == 0) or (len(uniprot_ids_b) == 0 and len(ncbi_ids_b) == 0):
        fail_to_find_uniprot_id += 1
        continue
    
    seq_a = get_seq(uniprot_ids_a, ncbi_ids_a, uniprot_id2seq, ncbi_id2seq)
    seq_b = get_seq(uniprot_ids_b, ncbi_ids_b, uniprot_id2seq, ncbi_id2seq)
    if seq_a is None or seq_b is None:
        fail_to_find_seqs_cnt += 1
        if seq_a is None:
            fail_to_find_seqs_uniprot.update(uniprot_ids_a)
            fail_to_find_seqs_ncbi.update(ncbi_ids_a)
            # print(uniprot_id_a, ncbi_id_a)
        if seq_b is None:
            fail_to_find_seqs_uniprot.update(uniprot_ids_b)
            fail_to_find_seqs_ncbi.update(ncbi_ids_b)
            # print(uniprot_id_b, ncbi_id_b)
        continue
        
    interaction_feats = {}
    for key in interaction_attr:
        interaction_feats[interaction_attr[key]] = row_dict[key]
        
    interactor_feats_a = {}
    interactor_feats_a["uniprot_ids"] = "|".join(uniprot_ids_a)
    interactor_feats_a["ncbi_ids"] = "|".join(ncbi_ids_a)
    for key in interactor_attr:
        interactor_feats_a[interactor_attr[key]] = row_dict[key + " A"]
        
    interactor_feats_b = {}
    interactor_feats_b["uniprot_ids"] = "|".join(uniprot_ids_b)
    interactor_feats_b["ncbi_ids"] = "|".join(ncbi_ids_b)
    for key in interactor_attr:
        interactor_feats_b[interactor_attr[key]] = row_dict[key + " B"]

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
print(f"fail_to_find_seqs: {fail_to_find_seqs_cnt}")
with open("fail_to_find_uniprot_id.log", "w") as fout:
    for id in fail_to_find_seqs_uniprot:
        fout.write(f"{id}\n")
with open("fail_to_find_ncbi_id.log", "w") as fout:
    for id in fail_to_find_seqs_ncbi:
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
