import json
import numpy as np
from collections import defaultdict, Counter
import random
from pathlib import Path

def output_statistics(dataset):
    print(f"Read {len(dataset)} data")
    antigen_set = set()
    vh_set = set()
    vl_set = set()
    pairs = set()
    pos_num = 0
    neg_num = 0
    
    for data in dataset:
        antigens = []
        vhs = []
        vls = []
        for interactor in data["interactors"]:
            if interactor["role"] == "antigen":
                antigens.append(interactor["seq"])
            elif interactor["role"] == "vh":
                vhs.append(interactor["seq"])
            elif interactor["role"] == "vl":
                vls.append(interactor["seq"])
        vh_set.update(vhs)
        vl_set.update(vls)
        antigen_set.update(antigens)
        if data["binding"]:
            pos_num += 1
        else:
            neg_num += 1
        vhs.sort()
        vls.sort()
        antigens.sort()
        pairs.add((tuple(vhs), tuple(vls), tuple(antigens)))
    print(f"Pos: {pos_num}; neg: {neg_num}; hit ratio {pos_num/len(dataset)}")
    print(f"Number of antigens: {len(antigen_set)}")
    print(f"Number of VHs: {len(vh_set)}")
    print(f"Number of VLs: {len(vl_set)}")
    print(f"Number of pairs: {len(pairs)}")
    print("Binding", Counter([x["binding"] for x in dataset]))

def split_by_binding(dataset, threshold=None):
    antigen2data = defaultdict(list)
    for data in dataset:
        binding = data["binding"]
        if threshold is not None:
            binding = binding > threshold
        antigen2data[binding].append(data)
    print(f"{len(antigen2data)} antigens")
    return antigen2data

def split_by_antigens(dataset):
    antigen2data = defaultdict(list)
    for data in dataset:
        antigens = []
        for interactor in data["interactors"]:
            if interactor["role"] == "antigen":
                antigens.append(interactor["id"])
        assert len(antigens) == 1
        antigen = antigens[0]
        antigen2data[antigen].append(data)
    print(f"{len(antigen2data)} antigens")
    # for antigen in antigen2data:
    #     print(antigen)
    #     output_statistics(antigen2data[antigen])
    return antigen2data

def random_subsample(dataset, subsample_num=650):
    rng = np.random.default_rng(seed=1005)
    # Sample 3 elements without replacement
    subsample_num = min(subsample_num, len(dataset))
    sampled_dataset = rng.choice(dataset, size=subsample_num, replace=False)
    return sampled_dataset


## Porebski HER2
split_by = "antigen"
binding_thres = 1 # 
data_path = Path("/data/rbg/users/wxsh/esm3/data/antibodies/porebski_et_al_2023_datasets/processed/porebski_her2_affmat/porebski_her2_affmat.json")

dataset = json.load(open(data_path))
output_statistics(dataset)
if split_by == "antigen":
    antigen2data = split_by_antigens(dataset)
elif split_by == "binding":
    antigen2data = split_by_binding(dataset, 0)
all_subsampled_dataset = []
for antigen in antigen2data:
    pos_data = [data for data in antigen2data[antigen] if data["binding"] > binding_thres]
    neg_data = [data for data in antigen2data[antigen] if data["binding"] <= binding_thres]
    print(f">>>{antigen}\nPos: {len(neg_data)}, neg: {len(neg_data)}")
    subsample_neg = random_subsample(neg_data, len(pos_data) * 9)
    all_subsampled_dataset.extend(pos_data)
    all_subsampled_dataset.extend(subsample_neg)
print(len(all_subsampled_dataset))
output_statistics(all_subsampled_dataset)

saving_data_dir = data_path.parent / "test_set_subsample_neg" / f"{data_path.stem}_subsample_{len(all_subsampled_dataset)}"
saving_data_dir.mkdir(exist_ok=True, parents=True)

json.dump(all_subsampled_dataset, open(saving_data_dir / f"data.json", "w"))

vh_id2seq = dict()
vl_id2seq = dict()
antigen_id2seq = dict()
for data in all_subsampled_dataset:
    for interactor in data["interactors"]:
        if interactor["role"] == "vh":
            vh_id2seq[interactor["id"]] = interactor["seq"]
        elif interactor["role"] == "vl":
            vl_id2seq[interactor["id"]] = interactor["seq"]
        elif interactor["role"] == "antigen":
            antigen_id2seq[interactor["id"]] = interactor["seq"]

def save_sequence_subsets(id2seq, path):
    with open(path, "w") as fout:
        for id, seq in id2seq.items():
            fout.write(f">{id}\n{seq}\n")

if len(vh_id2seq) > 0:
    save_sequence_subsets(vh_id2seq, saving_data_dir / f"vh_seqs.fasta")
if len(vl_id2seq) > 0:
    save_sequence_subsets(vl_id2seq, saving_data_dir / f"vl_seqs.fasta")
if len(antigen_id2seq) > 0:
    save_sequence_subsets(antigen_id2seq, saving_data_dir / f"antigen_seqs.fasta")