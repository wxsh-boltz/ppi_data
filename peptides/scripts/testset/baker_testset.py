import json
from collections import defaultdict, Counter
import json
from tqdm import tqdm
import numpy as np
from Bio import SeqIO

def read_alignment(path, identity_threshold=0.7, skip_self=False, swap_query_and_key=False, mismatch_threshold=None):
    query_id2similar_target_ids = defaultdict(list)
    with open(path) as fin:
        for line in fin:
            query_id, target_id, identity = line.strip().split()[:3]
            mismatch = int(line.strip().split()[4]) # 
            if swap_query_and_key:
                target_id, query_id = query_id, target_id
            identity = float(identity)
            if identity < identity_threshold:
                continue
            if mismatch_threshold is not None and mismatch > mismatch_threshold:
                continue
            if skip_self and query_id == target_id:
                continue
            query_id2similar_target_ids[query_id].append(target_id)
    return query_id2similar_target_ids

def entry_id2chain_ids(rcsb_ids):
    # rcsb_ids: List of string, like xxxx_A, xxxx_B
    r = defaultdict(list)
    for id in rcsb_ids:
        r[id.split("_")[0]].append(id.split("_")[1])
    return r


binder_si_threshold=0.7
target_si_threshold=0.7

# # Nature 2022
target_align_path_mode_2 = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/alnRes/target_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_2.m8"
target_align_path_mode_1 = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/alnRes/target_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_1.m8"
binder_align_path = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/alnRes/binder_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_1.m8"
binder_self_align_path = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/alnRes/binder_seqs_TO_binder_seqs.m8"
target_self_align_path = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/alnRes/target_seqs_TO_target_seqs.m8"
data_path = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/mini_proteins_cao_nature_2022.json"

# # Nature 2017
# target_align_path_mode_2 = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/alnRes/target_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_2.m8"
# target_align_path_mode_1 = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/alnRes/binder_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_1.m8"
# binder_align_path = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/alnRes/binder_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_1.m8"
# binder_self_align_path = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/alnRes/binder_seqs_TO_binder_seqs.m8"
# target_self_align_path = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/alnRes/target_seqs_TO_target_seqs.m8"
# data_path = "/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/mini_proteins_chevalier_nature_2017.json"

saving_path = data_path.replace(".json", f".test_target_{target_si_threshold}_binder_{binder_si_threshold}.json")


target_id2similar_rcsb_ids_mode_2 = read_alignment(target_align_path_mode_2, identity_threshold=target_si_threshold)
target_id2similar_rcsb_ids_mode_1 = read_alignment(target_align_path_mode_1, identity_threshold=target_si_threshold)
binder_id2similar_rcsb_ids = read_alignment(binder_align_path, identity_threshold=binder_si_threshold) # , mismatch_threshold=3
binder_id2binder_ids = read_alignment(binder_self_align_path, identity_threshold=1.0) # , mismatch_threshold=3
target_id2target_ids = read_alignment(target_self_align_path, identity_threshold=1.0)
print(f"len(binder_id2similar_rcsb_ids): {len(binder_id2similar_rcsb_ids)}")

all_data = json.load(open(data_path))
print(f"Read {len(all_data)} interactions")
pos_data = [data for data in all_data if data["binding"] > 0]
neg_data = [data for data in all_data if data["binding"] == 0]
print(f"pos_data: {len(pos_data)}")
print(f"neg_data: {len(neg_data)}")

data_not_in_pdb = []
all_targets = set()
testset_targets = set()
dataset_target2total_reads = dict()
binder_exclude_by_mmseqs = set()
target_exclude_by_mmseqs = set()
for i, data in enumerate(tqdm(all_data)):
    assert len(data["interactors"]) == 2
    for interactor in data["interactors"]:
        if interactor["role"] == "target":
            target = interactor
            all_targets.add(target["id"])
        elif interactor["role"] == "binder":
            binder = interactor
    
    if binder["id"] not in binder_id2binder_ids:
        binder_exclude_by_mmseqs.add((binder["id"], binder["seq"]))
        continue
    if target["id"] not in target_id2target_ids:
        target_exclude_by_mmseqs.add((target["id"], target["seq"]))
        continue

    binder_similar_rcsb_ids = entry_id2chain_ids(binder_id2similar_rcsb_ids[binder["id"]])
    target_similar_rcsb_ids_mode_2 = entry_id2chain_ids(target_id2similar_rcsb_ids_mode_2[target["id"]])
    target_similar_rcsb_ids_mode_1 = entry_id2chain_ids(target_id2similar_rcsb_ids_mode_1[target["id"]])
    overlap_entries_1 = set(target_similar_rcsb_ids_mode_2.keys()) & set(binder_similar_rcsb_ids.keys())
    overlap_entries_2 = set(target_similar_rcsb_ids_mode_1.keys()) & set(binder_similar_rcsb_ids.keys())
    overlap_entries = overlap_entries_1 | overlap_entries_2
    if len(overlap_entries) == 0:
        data_not_in_pdb.append(data)
    testset_targets.add(target["id"])

print(f"data_not_in_pdb: {len(data_not_in_pdb)}")
pos_data_not_in_pdb = [data for data in data_not_in_pdb if data["binding"] > 0]
neg_data_not_in_pdb = [data for data in data_not_in_pdb if data["binding"] == 0]
print(f"pos_data_not_in_pdb: {len(pos_data_not_in_pdb)} ({len(pos_data_not_in_pdb) / len(pos_data)})")
print(f"neg_data_not_in_pdb: {len(neg_data_not_in_pdb)} ({len(neg_data_not_in_pdb) / len(neg_data)})")
print(f"binder_exclude_by_mmseqs: {len(binder_exclude_by_mmseqs)}")
print(f"target_exclude_by_mmseqs: {len(target_exclude_by_mmseqs)}")
print(list(binder_exclude_by_mmseqs)[:3])
print(list(target_exclude_by_mmseqs)[:3])

print(f"all_targets: {len(all_targets)}")
print(f"testset_targets: {len(testset_targets)}: {testset_targets}")
new_targets = []
for target in all_targets:
    if len(set(target_id2similar_rcsb_ids_mode_1[target])) == 0 and len(set(target_id2similar_rcsb_ids_mode_2[target])) == 0:
        new_targets.append(target)
print("New targets:", len(new_targets), new_targets)

json.dump(data_not_in_pdb, open(saving_path, "w"))