from pathlib import Path
import json
from Bio import SeqIO

# simple concat
root_dir = Path("processed")
output_root_dir = Path("processed_merged")

output_root_dir.mkdir(exist_ok=True, parents=True)

dataset_name = "mini_proteins_cao_nature_2022"

# all_targets = []
target_id2seq = dict()
# all_binders = []
binder_id2seq = dict()
pairs = []

def get_id2seq(path, id2seq, id_prefix=None):
    for record in SeqIO.parse(path, "fasta"):
        if id_prefix is not None:
            new_id = f"{id_prefix}_{record.id}"
        else:
            new_id = record.id
        if new_id not in id2seq:
            id2seq[new_id] = str(record.seq)
        else:
            assert id2seq[new_id] == str(record.seq), f"Conflict ID: {new_id}, {path}"


for path in root_dir.glob("*.csv"):
    print(path)
    target_name = path.stem
    json_path = root_dir / f"{target_name}_exclude_confict_sc.json"
    # target_fasta_path = root_dir / f"{target_name}_target_unique_seqs.fasta"
    # binder_fasta_path = root_dir / f"{target_name}_binder_unique_seqs.fasta"
    dataset = json.load(open(json_path))

    for data in dataset:
        target, binder = None, None
        assert len(data["interactors"]) == 2
        for interactor in data["interactors"]:
            if interactor["role"] == "target":
                target = interactor
            elif interactor["role"] == "binder":
                binder = interactor
        assert target is not None and binder is not None

        new_binder_id = f"{target_name}_{binder['id']}"
        
        target_id2seq[target["id"]] = target["seq"]
        binder_id2seq[new_binder_id] = binder["seq"]
        binder["id"] = new_binder_id
        data["assay_id"] = f"{dataset_name}_{target_name}"
        
    # get_id2seq(target_fasta_path, target_id2seq)
    # get_id2seq(binder_fasta_path, binder_id2seq, target_name)
    # all_targets.extend(list([x for x in SeqIO.parse(target_fasta_path, "fasta")]))
    # all_binders.extend(list([x for x in SeqIO.parse(binder_fasta_path, "fasta")]))
    pairs.extend(dataset)

print(f"len(all_targets): {len(target_id2seq)}")
print(f"len(all_binders): {len(binder_id2seq)}")
print(f"len(pairs): {len(pairs)}")
print(f"Binder seqs {len(set(binder_id2seq.values()))}")

print(pairs[0])

json.dump(pairs, open(output_root_dir / f"{dataset_name}.json", "w"))
with open(output_root_dir / f"{dataset_name}_target_seqs.fasta", "w") as fout:
    for target_id in target_id2seq:
        fout.write(f">{target_id}\n{target_id2seq[target_id]}\n")
with open(output_root_dir / f"{dataset_name}_binder_seqs.fasta", "w") as fout:
    for binder_id in binder_id2seq:
        fout.write(f">{binder_id}\n{binder_id2seq[binder_id]}\n")
