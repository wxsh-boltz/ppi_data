
import pandas as pd
from collections import Counter, defaultdict
from pathlib import Path

root_dir = Path("/data/scratch/wxsh/data/peptides/supplemental_files/ngs_analysis/affinities")

sc_files = list(root_dir.glob("*.sc"))
skip_targets = ["IGF1R", "IL7Ra_graft", "H3_ssm_2"] # no designs provided 
# IL7Ra_graft: not sure
seq_pairs2designs = defaultdict(list)

# records = []
target_sequences = defaultdict(set)
binder_sequences = set()

for sc_file in sc_files:
    target = sc_file.stem

    if target.endswith("_ssm"):
        continue

    if target in skip_targets:
        continue

    path = f"/data/scratch/wxsh/data/peptides/supplemental_files/ngs_analysis/affinities/{target}.sc"
    # sequence_list_path = f"/data/scratch/wxsh/data/peptides/supplemental_files/sorting_ngs_data/{target}/sequences.list"
    # pdb_dir = Path(f"/data/rbg/shared/datasets/boltzgen_cao2022/supplemental_files/design_models_pdb/{target}")
    design_seq_path = f"/data/rbg/shared/datasets/boltzgen_cao2022/supplemental_files/design_models_sequence/{target}.seq"

    design_seqs = dict()
    with open(design_seq_path) as fin:
        for line in fin:
            try:
                peptide, protein, id = line.strip().split()
                design_seqs[id] = (peptide, protein)
            except Exception as e:
                print(line.strip())
                            
    print(target)

    df_binding = pd.read_csv(path, sep=" ")

    target_proteins = []
    binder_proteins = []
    for description in df_binding["description"]:
        peptide_protein, target_protein = design_seqs[description]
        target_proteins.append(target_protein)
        binder_proteins.append(peptide_protein)
    df_binding["target_protein"] = target_proteins
    df_binding["binder_protein"] = binder_proteins

    df_binding.to_csv(f"processed/{target}.csv",index=False)