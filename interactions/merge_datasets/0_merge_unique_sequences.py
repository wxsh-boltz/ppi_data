from Bio import SeqIO
import json
from collections import defaultdict
# fasta_paths = [
#     "/data/rbg/users/wxsh/esm3/data/interactions/intact/sequences_unique.fasta",
#     "/data/rbg/users/wxsh/esm3/data/interactions/bio_grid/process_scripts_boltz/sequences_unique.fasta"
# ]

max_seq_length = 4096

json_paths = [
    "/data/rbg/users/wxsh/esm3/data/interactions/bio_grid/BIOGRID-MV-Physical-4.4.245.direct_interactions.json",
    "/data/rbg/users/wxsh/esm3/data/interactions/intact/all_20250613/intact_direct_interactions.json"
]

saving_path = f"/data/rbg/users/wxsh/esm3/data/interactions/intact_and_biogrid_direct_interactions/sequences_unique_maxL_{max_seq_length}.fasta"

seq2ids = defaultdict(list)

for path in json_paths:
    records = json.load(open(path))
    print(f"Read {len(records)} interactions from {path}.")
    for record in records:
        interactors = record["interactors"]
        for interactor in interactors:
            if len(interactor["seq"]) <= max_seq_length:
                seq2ids[interactor["seq"]].append(interactor["id"].replace("intact:", "intact_"))

print(f"Read {len(seq2ids)} sequences.")

print(f"Max sequence length: {max([len(x) for x in seq2ids])}")
# print(f"Max sequence length: {max([len(x) for in seq2ids])}")


with open(saving_path, "w") as fout:
    for seq in seq2ids:
        ids = set(seq2ids[seq])
        ids = list(ids)
        ids.sort()
        if len(ids) > 1:
            id_others = "|".join(ids[1:])
        else:
            id_others = ""
        # {id_others}
        fout.write(f">{ids[0]}\n{seq}\n")