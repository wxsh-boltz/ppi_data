from Bio import SeqIO
import json
from collections import defaultdict

json_paths = {
    "biogrid": "/data/rbg/users/wxsh/esm3/data/interactions/bio_grid/BIOGRID-MV-Physical-4.4.245.direct_interactions.json",
    "intact": "/data/rbg/users/wxsh/esm3/data/interactions/intact/all_20250613/intact_direct_interactions.json"
}

merged_seq_path = f"/data/rbg/users/wxsh/esm3/data/interactions/intact_and_biogrid_direct_interactions/sequences_unique_maxL_4096.fasta"
saving_path = f"/data/rbg/users/wxsh/esm3/data/interactions/intact_and_biogrid_direct_interactions/pos_maxL_4096.json"

seq2new_id = dict()
for record in SeqIO.parse(merged_seq_path, "fasta"):
    seq2new_id[str(record.seq)] = record.id

print(f"Number of unique seqs: {len(seq2new_id)}")

seq_pair2attributes = defaultdict(list)
for dataset_name, path in json_paths.items():
    records = json.load(open(path))
    print(f"Read {len(records)} interactions from {path} for {dataset_name}.")
    for record in records:
        interactors = record["interactors"]
        assert len(interactors) == 2
        assert interactors[0]["seq"] <= interactors[1]["seq"]
        # interactors = sorted(interactors, key=lambda x: x["seq"])
        for attr in record["attributes"]:
            attr["dataset_name"] = dataset_name
        seq_pair2attributes[tuple([x["seq"] for x in interactors])].extend(record["attributes"])

print(f"Merged data size: {len(seq_pair2attributes)}")

new_records = []
exclude_number_pairs_for_too_long_seqs = 0
for (seq1, seq2), attributes in seq_pair2attributes.items():
    seq1_id = seq2new_id.get(seq1)
    seq2_id = seq2new_id.get(seq2)
    if seq1_id is None or seq2_id is None: # exclude because it is too long
        # print(len(seq1), len(seq2))
        exclude_number_pairs_for_too_long_seqs += 1
        continue
    new_interactors = [{"seq": seq1, "id": seq1_id}, {"seq": seq2, "id": seq2_id}]
    new_records.append({"interactors": new_interactors, "attributes": attributes, "binding": True})

print(f"Merged data size (after filtering): {len(new_records)}. Exclude {exclude_number_pairs_for_too_long_seqs} for too long seqs.")

json.dump(new_records, open(saving_path, "w"))