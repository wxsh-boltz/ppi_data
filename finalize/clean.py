import json
import numpy as np
from collections import Counter

input_path = "/data/rbg/users/wxsh/esm3/data/peptides/ppikd/processed/data_train/data.pos_and_neg.json"
input_path = "/data/rbg/users/wxsh/esm3/data/peptides/PANCS_Binders/processed/pancs_merged_train/data_train/data_subsample_17511/data.json"
input_path = "/data/rbg/users/wxsh/esm3/data/interactions/intact_and_biogrid_direct_interactions/pos_maxLL_1500.plddt70/data_train/data.pos_and_neg.json"
input_path = "/data/rbg/users/wxsh/esm3/data/skempi/processed/skempi_v2_by_pdb.train_train/data.json"

save_path = input_path.replace(".json", ".clean.json")
print(input_path)
dataset = json.load(open(input_path))
print(len(dataset))

def quality_control_seq(seq, max_ratio_of_x = 0.5):
    seq = seq.upper()
    if np.mean(np.asarray(list(seq)) == "X") >= max_ratio_of_x:
        return None
    return seq.upper()

clean_dataset = []
for data in dataset:
    _skip = False
    new_interactors = []
    for interactor in data["interactors"]:
        seq_cananical = quality_control_seq(interactor["seq"])
        if seq_cananical is None:
            _skip=True
            break
        else:
            interactor["seq"] = seq_cananical
    if _skip:
        continue
    clean_dataset.append(data)
print(len(clean_dataset))

for data in clean_dataset:
    for interactor in data["interactors"]:
        assert any(ch.isupper() for ch in interactor["seq"])

print(Counter([data["binding"] for data in clean_dataset]).most_common())
json.dump(clean_dataset, open(save_path, "w"))