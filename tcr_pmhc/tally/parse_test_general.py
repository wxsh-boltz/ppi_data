from pathlib import Path
import pandas as pd
from collections import Counter, defaultdict
import json
import numpy as np

saving_parent_dir = Path("processed")
saving_parent_dir.mkdir(exist_ok=True, parents=True)

dataset_name2path = {
    # "jacks": "/data/rbg/shared/projects/boltz_tcr/datasets/Tyler/parsed_data/filtered_manifest_no_crossreactive20250223.csv",
    "jacks_v2": "/data/rbg/shared/projects/boltz_tcr/datasets/Tyler/parsed_data/downsampled_manifest_20250511.csv",
    "immrep23": "/data/rbg/shared/projects/boltz_tcr/datasets/immrep2023/processed/test_restitched.csv",
    "immrep25_train_no10x": "/data/rbg/shared/projects/boltz_tcr/datasets/immrep25/processed/balanced_train_20250315.csv",
    "immrep25_val_no10xv2": "/data/rbg/shared/projects/boltz_tcr/datasets/immrep25/processed/balanced_val_20250315.csv",
    "mouse_mixed_val": "/data/rbg/shared/projects/boltz_tcr/datasets/vdjdb/processed/finetune/balanced_val_20250327_mouse_mh2d1_h2kb.csv",
    "mouse_h2kb_val": "/data/rbg/shared/projects/boltz_tcr/datasets/vdjdb/processed/finetune/balanced_val_20250327_mouse_h2kb.csv",
    "mouse_mh2d1_val": "/data/rbg/shared/projects/boltz_tcr/datasets/vdjdb/processed/finetune/balanced_val_20250327_mouse_mh2d1.csv",
    "mouse_mixed_train": "/data/rbg/shared/projects/boltz_tcr/datasets/vdjdb/processed/finetune/balanced_train_20250327_mouse_mh2d1_h2kb.csv",
    "mouse_h2kb_train": "/data/rbg/shared/projects/boltz_tcr/datasets/vdjdb/processed/finetune/balanced_train_20250327_mouse_h2kb.csv",
    "mouse_mh2d1_train": "/data/rbg/shared/projects/boltz_tcr/datasets/vdjdb/processed/finetune/balanced_train_20250327_mouse_mh2d1.csv",
}

dataset_name2cols = {
    "immrep23": ["Peptide", "va", "vb", "mhc.a.seq.cropped", "Label"],
    # "jacks": ["peptide", "TCR_VA", "TCR_VB", "MHC_A", "bind"],
    "jacks_v2": ["peptide", "TCR_VA", "TCR_VB", "MHC_A", "bind"],
    "immrep25_train_no10x": ["Peptide", "TCRa.seq", "TCRb.seq", "mhc.a.seq.cropped", "bind"],
    "immrep25_val_no10xv2": ["Peptide", "TCRa.seq", "TCRb.seq", "mhc.a.seq.cropped", "bind"],
    "mouse_mixed_val": ["antigen.epitope", "va", "vb", "mhc.a.seq.cropped", "bind"],
    "mouse_h2kb_val": ["antigen.epitope", "va", "vb", "mhc.a.seq.cropped", "bind"],
    "mouse_mh2d1_val": ["antigen.epitope", "va", "vb", "mhc.a.seq.cropped", "bind"],
    "mouse_mixed_train": ["antigen.epitope", "va", "vb", "mhc.a.seq.cropped", "bind"],
    "mouse_h2kb_train": ["antigen.epitope", "va", "vb", "mhc.a.seq.cropped", "bind"],
    "mouse_mh2d1_train": ["antigen.epitope", "va", "vb", "mhc.a.seq.cropped", "bind"],
}

dataset_name2positive_only = {
    "immrep23": False,
    # "jacks": False,
    "jacks_v2": False,
    "immrep25_train_no10x": True,
    "immrep25_val_no10xv2": True,
    "mouse_mixed_val": True,
    "mouse_h2kb_val": True,
    "mouse_mh2d1_val": True,
    "mouse_mixed_train": True,
    "mouse_h2kb_train": True,
    "mouse_mh2d1_train": True,

}

solutions_df = pd.read_csv(
    "/data/rbg/shared/projects/boltz_tcr/datasets/immrep2023/raw/solutions.csv"
)

def check_seq_lengths(seqs):
    lengths = np.asarray([len(s) for s in seqs])
    print(f"Average L: {np.mean(lengths)}, max L: {np.max(lengths)}, min L: {np.min(lengths)}")
    
for dataset_name, path in dataset_name2path.items():
    print(f">>> Processing {dataset_name} from {path}.")
    peptide_col, va_col, vb_col, mhc_col, label_col = dataset_name2cols[dataset_name]

    df = pd.read_csv(path) # .reset_index()
    print(f"Read {len(df)} data from {dataset_name}.")
    # df.set_index("name", inplace=True)
    if dataset_name in ("immrep23",):
        df = df.merge(
            solutions_df[["ID", "Label"]], left_on="ID", right_on="ID", how="left"
        )
    
    df = df.dropna(subset=[peptide_col, va_col, vb_col, mhc_col, label_col])
    
    print(f"Read {len(df)} data from {dataset_name} (after dropna).")
    print("Label:", Counter(df[label_col]))
    print("Peptide:", Counter(df[peptide_col]).most_common(3))

    mhc_all = set()
    tcr_va_all = set()
    tcr_vb_all = set()
    peptide_all = set()

    pairs = defaultdict(list)
    for row in df.to_dict(orient='records'):
        mhc = row.pop(mhc_col)
        tcr_va = row.pop(va_col)
        tcr_vb = row.pop(vb_col)
        peptide = row.pop(peptide_col)
        mhc_all.add(mhc)
        tcr_va_all.add(tcr_va)
        tcr_vb_all.add(tcr_vb)
        peptide_all.add(peptide)
        pairs[(mhc, tcr_va, tcr_vb, peptide)].append(row)
    print(f"Pairs: {len(pairs)}")
    print(f"MHC_A: {len(mhc_all)}")
    check_seq_lengths(mhc_all)
    print(f"TCR_VA: {len(tcr_va_all)}")
    check_seq_lengths(tcr_va_all)
    print(f"TCR_VB: {len(tcr_vb_all)}")
    check_seq_lengths(tcr_vb_all)
    print(f"peptide: {len(peptide_all)}")
    check_seq_lengths(peptide_all)

    # dataset2df[dataset_name] = df
    saving_dir = saving_parent_dir / dataset_name
    saving_dir.mkdir(exist_ok=True, parents=True)

    mhc_all = list(mhc_all)
    mhc_all.sort()
    mhc_seq2id = dict()
    with open(saving_dir / "mhc_seqs.fasta", "w") as fout:
        for i, seq in enumerate(mhc_all):
            mhc_seq2id[seq] = f"{dataset_name}_mhc_{i}"
            fout.write(f">{mhc_seq2id[seq]}\n{seq}\n")

    tcr_va_all = list(tcr_va_all)
    tcr_va_all.sort()
    tcr_va_seq2id = dict()
    with open(saving_dir / "tcr_va_seqs.fasta", "w") as fout:
        for i, seq in enumerate(tcr_va_all):
            tcr_va_seq2id[seq] = f"{dataset_name}_tcr_va_{i}"
            fout.write(f">{tcr_va_seq2id[seq]}\n{seq}\n")

    tcr_vb_all = list(tcr_vb_all)
    tcr_vb_all.sort()
    tcr_vb_seq2id = dict()
    with open(saving_dir / "tcr_vb_seqs.fasta", "w") as fout:
        for i, seq in enumerate(tcr_vb_all):
            tcr_vb_seq2id[seq] = f"{dataset_name}_tcr_vb_{i}"
            fout.write(f">{tcr_vb_seq2id[seq]}\n{seq}\n")

    peptide_all = list(peptide_all)
    peptide_all.sort()
    peptide_seq2id = dict()
    with open(saving_dir / "peptide_seqs.fasta", "w") as fout:
        for i, seq in enumerate(peptide_all):
            peptide_seq2id[seq] = f"{dataset_name}_peptide_{i}"
            fout.write(f">{peptide_seq2id[seq]}\n{seq}\n")

    output_list = []
    for (mhc, tcr_va, tcr_vb, peptide), attrs in pairs.items():
        interactors = []
        interactors.append({"seq": peptide, "id": peptide_seq2id[peptide], "role": "peptide"})
        interactors.append({"seq": tcr_va, "id": tcr_va_seq2id[tcr_va], "role": "tcr_va"})
        interactors.append({"seq": tcr_vb, "id": tcr_vb_seq2id[tcr_vb], "role": "tcr_vb"})
        interactors.append({"seq": mhc, "id": mhc_seq2id[mhc], "role": "mhc"})
        binding = set()
        for attribute in attrs:
            binding.add(attribute[label_col])
        if len(binding) > 1:
            print(binding)
            continue
        binding = list(binding)[0]
        if isinstance(binding, str):
            if binding == "bind":
                binding_bool = True
            else:
                binding_bool = False
        elif isinstance(binding, float) or isinstance(binding, int):
            binding_bool = bool(binding)
        else:
            raise ValueError(type(binding))

        output_list.append({
            "interactors": interactors, 
            "binding": binding_bool, 
            "attributes": attrs,
            "assay_id": dataset_name}
            )
        
    json.dump(output_list, open(saving_dir / f"{dataset_name}.json", "w"))

    if dataset_name2positive_only[dataset_name]:
        output_list_pos = [d for d in output_list if d["binding"]]
        json.dump(output_list_pos, open(saving_dir / f"{dataset_name}.pos.json", "w"))

    print("Done.\n")