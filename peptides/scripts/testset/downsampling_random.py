import json
import numpy as np
from collections import defaultdict, Counter
import random
from pathlib import Path

split_by = "target" 
BIND_THRES = 0 # > BIND_THRES: binder, <= BIND_THRES: non-binders
MAX_BINDERS = 100000000000
data_path = Path("/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/mini_proteins_cao_nature_2022.test_target_0.7_binder_0.7.json")
saving_data_dir = Path("/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/test_set_subsample")

# BIND_THRES = 2 # > BIND_THRES: binder, <= BIND_THRES: non-binders
# MAX_BINDERS = 300
# data_path = Path("/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/mini_proteins_chevalier_nature_2017.test_target_0.7_binder_0.7.json")
# saving_data_dir = Path("/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/testset_subsample")

### Pancs
BIND_THRES = 0 # > BIND_THRES: binder, <= BIND_THRES: non-binders
MAX_BINDERS = 100000000000
data_path = Path("/data/rbg/users/wxsh/esm3/data/peptides/PANCS_Binders/processed/pancs_merged.test.json")
saving_data_dir = Path("/data/rbg/users/wxsh/esm3/data/peptides/PANCS_Binders/processed/testset_subsample")
split_by = "assay" # or assay

def output_statistics(dataset):
    print(f"Read {len(dataset)} data")
    target_set = set()
    binder_set = set()
    pairs = set()
    pos_num = 0
    neg_num = 0
    print(f"Binding: {Counter([data['binding'] for data in dataset])}")
    
    for data in dataset:
        targets = []
        binders = []
        for interactor in data["interactors"]:
            if interactor["role"] == "target":
                targets.append(interactor["seq"])
            elif interactor["role"] == "binder":
                binders.append(interactor["seq"])
        if data["binding"] > BIND_THRES:
            pos_num += 1
        else:
            neg_num += 1
        assert len(targets) == 1
        assert len(binders) == 1
        target_set.add(targets[0])
        binder_set.add(binders[0])
        pairs.add((targets[0], binders[0]))
    print(f"Pos: {pos_num}; neg: {neg_num}; hit ratio {pos_num/len(dataset)}")
    print(f"Number of targets: {len(target_set)}")
    print(f"Number of binders: {len(binder_set)}")
    print(f"Number of pairs: {len(pairs)}")

def split_by_targets(dataset):
    target2data = defaultdict(list)
    for data in dataset:
        targets = []
        for interactor in data["interactors"]:
            if interactor["role"] == "target":
                targets.append(interactor["id"])
        assert len(targets) == 1
        target = targets[0]
        target2data[target].append(data)
    print(f"{len(target2data)} targets")
    return target2data

def split_by_assay(dataset):
    assay2data = defaultdict(list)
    for data in dataset:
        assay2data[data["assay_id"]].append(data)
    print(f"{len(assay2data)} assay")
    return assay2data

def random_subsample(dataset, subsample_num=650):
    rng = np.random.default_rng(seed=1005)
    subsample_num = min(subsample_num, len(dataset))
    sampled_dataset = rng.choice(dataset, size=subsample_num, replace=False)
    return sampled_dataset

saving_data_dir.mkdir(exist_ok=True, parents=True)

dataset = json.load(open(data_path))
output_statistics(dataset)
if split_by == "target":
    target2data = split_by_targets(dataset)
elif split_by == "assay":
    target2data = split_by_assay(dataset)

all_subsampled_dataset = []
for target in target2data:
    binders = [d for d in target2data[target] if d["binding"] > BIND_THRES]
    non_binders = [d for d in target2data[target] if d["binding"] <= BIND_THRES]
    
    print(f"Target/assay: {target}, binders: {len(binders)}, non_binders: {len(non_binders)}")
    if len(binders) > MAX_BINDERS:
        binders = random_subsample(binders, MAX_BINDERS)
    
    neg_subsample_num = len(binders) * 9
    subsample_non_binders = random_subsample(non_binders, neg_subsample_num)
    print(f"subsample_non_binders: {len(subsample_non_binders)}")
    all_subsampled_dataset.extend(binders)
    all_subsampled_dataset.extend(subsample_non_binders)
    
print(len(all_subsampled_dataset))
output_statistics(all_subsampled_dataset)

json.dump(all_subsampled_dataset, open(saving_data_dir / f"{data_path.stem}_subsample_{len(all_subsampled_dataset)}.json", "w"))

target_id2seq = dict()
binder_id2seq = dict()
for data in all_subsampled_dataset:
    for interactor in data["interactors"]:
        if interactor["role"] == "target":
            target_id2seq[interactor["id"]] = interactor["seq"]
        elif interactor["role"] == "binder":
            binder_id2seq[interactor["id"]] = interactor["seq"]

def save_sequence_subsets(id2seq, path):
    with open(path, "w") as fout:
        for id, seq in id2seq.items():
            fout.write(f">{id}\n{seq}\n")

if len(binder_id2seq) > 0:
    save_sequence_subsets(binder_id2seq, saving_data_dir / f"{data_path.stem}_subsample_{len(all_subsampled_dataset)}_binder_seqs.fasta")
if len(target_id2seq) > 0:
    save_sequence_subsets(target_id2seq, saving_data_dir / f"{data_path.stem}_subsample_{len(all_subsampled_dataset)}_target_seqs.fasta")