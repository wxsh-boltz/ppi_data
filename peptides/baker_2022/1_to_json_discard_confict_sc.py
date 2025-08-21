import pandas as pd
from collections import defaultdict
import json
from pathlib import Path

root_dir = Path("/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed/")

for file in root_dir.glob("*.csv"):
    # print(file)
    target_name = file.stem #  "EGFR"
    print(target_name)

    path = f"/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed/{target_name}.csv"

    target_fasta_saving_path = path.replace(".csv", "_target_unique_seqs.fasta")
    binder_fasta_saving_path = path.replace(".csv", "_binder_unique_seqs.fasta")

    output_path = path.replace(".csv", "_exclude_confict_sc.json")

    df = pd.read_csv(path)

    seq_pairs2experiments = defaultdict(list)
    target_protein_seqs = set()
    binder_protein2ids = defaultdict(set)
    for idx, row in df.iterrows():
        row_dict = row.to_dict()
        target_protein = row_dict.pop("target_protein")
        binder_protein = row_dict.pop("binder_protein")
        target_protein_seqs.add(target_protein)
        binder_protein2ids[binder_protein].add(row_dict["description"])
        seq_pairs2experiments[(target_protein, binder_protein)].append(row_dict)


    target_protein_seqs = list(target_protein_seqs)
    target_protein_seqs.sort()
    target_protein_seq2unique_id = dict()
    with open(target_fasta_saving_path, "w") as fout:
        for i, seq in enumerate(target_protein_seqs):
            id_unique = f"{target_name}_{i}"
            target_protein_seq2unique_id[seq] = id_unique
            fout.write(f">{id_unique}\n{seq}\n")
    binder_protein_seq2unique_id = dict()
    with open(binder_fasta_saving_path, "w") as fout:
        for seq in binder_protein2ids:
            ids = list(binder_protein2ids[seq])
            ids.sort()
            # if len(ids) >= 2:
            #     print(ids)
            id_unique = ids[0]
            binder_protein_seq2unique_id[seq] = id_unique
            # fout.write(f">{id_unique} {'|'.join(ids[1:])}\n{seq}\n")
            fout.write(f">{id_unique}\n{seq}\n")

    ### 
    print(f"DF size: {len(df)}, seq pairs: {len(seq_pairs2experiments)}")
    print(f"Target seqs: {len(target_protein_seq2unique_id)}")
    print(f"Binder seqs: {len(binder_protein_seq2unique_id)}")

    records = []
    # records_exclude_conflicts = []

    for target_protein, binder_protein in seq_pairs2experiments:
        interactor = [
            {"id": target_protein_seq2unique_id[target_protein], "seq": target_protein, "role": "target"},
            {"id": binder_protein_seq2unique_id[binder_protein], "seq": binder_protein, "role": "binder"},
        ]
        # attributes = []
        sc_low = []
        sc_high = []
        for experiment in seq_pairs2experiments[(target_protein, binder_protein)]:
            # print(experiment)
            if "binder_4000_nm" in experiment:
                sc_high.append(experiment["binder_4000_nm"])
            
            if "binder_400_nm" in experiment:
                sc_low.append(experiment["binder_400_nm"])
            elif "binder_800_nm" in experiment: # For H3
                sc_low.append(experiment["binder_800_nm"])
        
        sc_low = set(sc_low)
        sc_high = set(sc_high)
        assert len(sc_low) <= 1 and len(sc_high) <= 1
        assert len(sc_low) == 1 or len(sc_high) == 1

        if len(sc_low) == 0:
            sc_low = False
        else:
            sc_low = list(sc_low)[0]
        if len(sc_high) == 0:
            sc_high = False
        else:
            sc_high = list(sc_high)[0]
        
        if sc_low and sc_high:
            binding = 2 # "strong"
        elif (not sc_low) and (not sc_high):
            binding = 0 # "no"
        elif (not sc_low) and (sc_high):
            binding = 1 # "weak"
        else:
            continue # skip this
        
        records.append({
            "interactors": interactor,
            "attributes": seq_pairs2experiments[(target_protein, binder_protein)],
            "binding": binding
            })
        # else:
            # binding = "no"


    print(f"Records: {len(records)}")
    json.dump(records, open(output_path, "w"))
