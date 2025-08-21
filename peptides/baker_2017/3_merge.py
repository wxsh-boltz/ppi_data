import json
from pathlib import Path
botn_path = "processed/botn.json"
ha_path = "processed/ha.json"
dataset_name = "mini_proteins_chevalier_nature_2017"

merge_dir = Path("processed_merged")
merge_dir.mkdir(exist_ok=True, parents=True)

botn = json.load(open(botn_path))
ha = json.load(open(ha_path))
print(f"botn: {len(botn)}, ha: {len(ha)}")

target_id2seq = {}
binder_id2seq = {}
merge_data = botn + ha
for data in merge_data:
    assert len(data["interactors"]) == 2
    target, binder = None, None
    for interactor in data["interactors"]:
        if interactor["role"] == "target":
            target = interactor
        elif interactor["role"] == "binder":
            binder = interactor
    assert target is not None
    assert binder is not None

    if target["id"] in target_id2seq:
        assert target_id2seq[target["id"]] == target["seq"], f"{target['id']}, {target_id2seq[target['id']]}, {target['seq']}"
    target_id2seq[target["id"]] = target["seq"]
    if binder["id"] in binder_id2seq:
        assert binder_id2seq[binder["id"]] == binder["seq"]
    binder_id2seq[binder["id"]] = binder["seq"]
    data["assay_id"] = f"{dataset_name}_{target['id'].split('_')[0]}"
print(len(target_id2seq), len(binder_id2seq))
print(len(set(binder_id2seq.values())))

print(f"Merged: {len(merge_data)}")
json.dump(merge_data, open(merge_dir/ f"{dataset_name}.json", "w"))

with open(merge_dir / f"{dataset_name}_target_seqs.fasta", "w") as fout:
    for target_id in target_id2seq:
        fout.write(f">{target_id}\n{target_id2seq[target_id]}\n")
with open(merge_dir / f"{dataset_name}_binder_seqs.fasta", "w") as fout:
    for binder_id in binder_id2seq:
        fout.write(f">{binder_id}\n{binder_id2seq[binder_id]}\n")