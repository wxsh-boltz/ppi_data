import json
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict, Counter
from pathlib import Path

all_seqs_path = "/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/pdb_05_08_2025_wx/sequences.json"
date_cutoff = np.datetime64("2023-06-01")
saving_dir = "polymers_after_20230601"
max_num_interactor = 2
dataset_name = "rcsb_2prot_after_20230601"
# manifest_path = "/data/rbg/shared/projects/foldeverything/boltz2/rcsb/targets/manifest.json"
# saving_path = "rcsb_boltz2_train_seqs.fasta"

dataset = json.load(open(all_seqs_path))
print(f"Read {len(dataset)} data")

complexes = defaultdict(list)
complexes_entry_count = 0
for entry in dataset:
    if np.datetime64(entry["released"]) > date_cutoff: # skip 
        if len(entry["proteins"]) > 1:
            if len(entry["rnas"]) > 0 or len(entry["dnas"]) > 0: # or len(entry["ligands"]) > 0
                continue
            protein_seqs = [protein[4]["seq_boltz_extended_and_gemmi"] for protein in entry["proteins"]]
            # protein_seqs_pdbx = [protein[4]["seq_pdbx"] for protein in entry["proteins"]]
            # print(protein_seqs)
            # print(protein_seqs_pdbx)
            if len(protein_seqs) > max_num_interactor:
                continue
            protein_seqs.sort()
            complexes[tuple(protein_seqs)].append(entry)
            complexes_entry_count += 1
            
print(len(complexes), complexes_entry_count)
# number_of_interactors = [len(seqs) for seqs in complexes]
# print(Counter(number_of_interactors).most_common())
# number_of_interactors = [len(set(list(seqs))) for seqs in complexes]
# print(Counter(number_of_interactors).most_common())
# for seqs in complexes:
#     if len(seqs) >= 5:
#         print(complexes[seqs])
#         exit()

records = []
id2seq = dict()
seq2id = dict()
# max_interactors = 10

for proteins_seqs in complexes:
    repr_entry = complexes[proteins_seqs][0]
    interactors = []
    for protein in repr_entry["proteins"]:
        pdb_id, chain_id, _, _, seq = protein
        seq = seq["seq_boltz_extended_and_gemmi"]
        seq_id = f"{pdb_id}_{chain_id}"
        assert seq_id not in id2seq, protein
        id2seq[seq_id] = seq
        if seq not in seq2id:
            seq2id[seq] = seq_id
        interactors.append({"id": seq2id[seq], "seq": seq})
    interactors = sorted(interactors, key=lambda x: x["seq"])
    
    # for protein_seq in proteins_seqs:
    all_pdb_ids = []
    for entry in complexes[proteins_seqs]:
        chain_seq2chain_id = {protein[4]["seq_boltz_extended_and_gemmi"]: f"{protein[0]}_{protein[1]}" for protein in entry["proteins"]}
        all_pdb_ids.append([chain_seq2chain_id[seq] for seq in proteins_seqs])

    # all_pdb_ids = [ for entry in complexes[proteins_seqs]]
    records.append({
        "interactors": interactors, 
        "attributes": {"pdb_ids": all_pdb_ids}, 
        "assay_id": dataset_name})

print(f"Read {len(id2seq)} seqs, but the unique seqs are only {len(set(id2seq.values()))}")
print(f"Number of records: {len(records)}")
print(f"Number of unique seqs: {len(seq2id)}")
print(records[0])
print(records[-1])



json.dump(records, open(Path(saving_dir) / f"{dataset_name}.json", "w"))

with open(Path(saving_dir) / f"{dataset_name}.fasta", "w") as fout:
    for seq, id in seq2id.items():
        fout.write(f">{id}\n{seq}\n")