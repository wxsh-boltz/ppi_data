import pandas as pd
from collections import Counter, defaultdict
from Bio import SeqIO
import json
from pathlib import Path
from Bio.PDB import PDBParser, PPBuilder
from tqdm import tqdm

import gemmi
from pathlib import Path
from tqdm import tqdm

def extract_protein_sequences(pdb_file):
    """
    Extract protein sequences for each chain from a PDB file using Gemmi.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        dict: Dictionary with chain IDs as keys and sequences as values
    """
    # Read the PDB file
    structure = gemmi.read_structure(pdb_file)
    structure.setup_entities()
    chain2seq = dict()
    for model in structure:
        for chain in model:
            # if not chain.is_polymer():
            #     continue

            residues = chain.get_polymer()
            seq = ""
            for res in residues:
                try:
                    aa = gemmi.find_tabulated_residue(res.name).one_letter_code
                    if not aa:
                        aa = "X"
                except Exception:
                    aa = "X"
                seq += aa

            chain2seq[chain.name] = seq
    return chain2seq


root_dir=Path("/data/scratch/wxsh/data/peptides/nature_2017_baker/archive_MassivelyParallelDeNovoProteinDesignForTargetedTherapeutics_2017/")
path = root_dir / "01_HighThroughputDesign_BotNTB_relax.dat"
pdb_dir = root_dir / Path("02_HighThroughputDesign_BotNTB_relax_models")
csv_path = Path("botn.csv")
dataset_name = "botn"

def parse(file):
    res = dict()
    
    for line in open(file):
        # print(line)
        if line.startswith("##"):
            continue
        if line.strip().startswith("#ID"):
            head_line = line.strip().split()
            # print(head_line)
            for col in head_line:
                res[col] = list()
            continue
        if line.startswith("#"):
            continue

        if line.strip():
            line = line.strip().split()
            for key, value in zip(head_line, line):
                res[key].append(value)
    return pd.DataFrame(res)

# def extract_protein_sequences(pdb_file):
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure("protein", pdb_file)
    
#     ppb = PPBuilder()
#     sequences = {}
    
#     for model in structure:
#         for chain in model:
#             peptides = ppb.build_peptides(chain)
#             if peptides:
#                 # Concatenate segments (in case of gaps)
#                 sequence = ''.join([str(pp.get_sequence()) for pp in peptides])
#                 sequences[chain.id] = sequence
#     return sequences

if not csv_path.exists():
    df_botn = parse(path)
    print(f"Len df_botn: {len(df_botn)}")
    pdb_files = list(pdb_dir.rglob("*.pdb"))
    print("pdb_files:", len(pdb_files))
    id2seq_pairs = dict()
    for pdb_file in tqdm(pdb_files):
        id = pdb_file.stem
        sequences = extract_protein_sequences(str(pdb_file))
        assert len(sequences) == 2, sequences
        if id in id2seq_pairs:
            assert id2seq_pairs[id] == (sequences["A"], sequences["B"])
        id2seq_pairs[id] = (sequences["A"], sequences["B"])

    target_seqs = []
    sequence_from_pdb = []
    for seq, cat, id in zip(df_botn["Sequence"], df_botn["Category"], df_botn["#ID"]):
        target, peptide = id2seq_pairs[id]
        target_seqs.append(target)
        sequence_from_pdb.append(peptide)
        
    df_botn["target_seq"] = target_seqs
    df_botn["sequence_from_pdb"] = sequence_from_pdb
    df_botn.to_csv(csv_path, index=False)
else:
    df_botn = pd.read_csv(csv_path)

print(f"df_botn: {len(df_botn)}")

# print(f"[Number of unique sequences]: peptide {len(set(df_botn['Sequence']))}, target {len(set(df_botn['target_seq']))}, peptide from pdb {len(set(df_botn['sequence_from_pdb']))}")

# seq2categories = defaultdict(list)
# peptide_seq2id = dict()
# all_target_seqs = set()
# for seq, cat, id, target, sequence_from_pdb in zip(df_botn["Sequence"], df_botn["Category"], df_botn["#ID"], df_botn["target_seq"], df_botn["sequence_from_pdb"]):
#     if seq not in peptide_seq2id:
#         peptide_seq2id[seq] = id
#     seq2categories[(seq, target)].append((cat, id, sequence_from_pdb))
#     # assert seq == sequence_from_pdb, f"{seq} {sequence_from_pdb}, {sequence_from_pdb in seq}"
#     all_target_seqs.add(target)
# print(f"[Number of seq pairs]: {len(seq2categories)}")

# count = dict(Counter(df_botn["Category"])) # .most_common()
# print("[Categories]", count[0], count[1], count[2], count[3], sum(count.values()))

# all_target_seqs = list(all_target_seqs)
# all_target_seqs.sort()
# target_seq2id = dict()
# with open(f"{dataset_name}_target_seqs.fasta", "w") as fout:
#     for i, seq in enumerate(all_target_seqs):
#         fout.write(f">{dataset_name}_target_{i}\n{seq}\n")
#         target_seq2id[seq] = f"{dataset_name}_target_{i}"
# with open(f"{dataset_name}_binder_seqs.fasta", "w") as fout:
#     for seq, id in peptide_seq2id.items():
#         fout.write(f">{id}\n{seq}\n")


# botn_records = []
# conflict_cases = 0
# for seq, target in seq2categories:
#     peptide_id = peptide_seq2id[seq]
    
#     attrs = []
#     categories = []
#     for cat, id, sequence_from_pdb in seq2categories[(seq, target)]:
#         categories.append(cat)
#         attrs.append({"category": cat, "id": id, "sequence_from_pdb": sequence_from_pdb})

#     categories_set = set(categories)
#     if len(categories_set) > 1:
#         print(categories_set)
#         conflict_cases += 1
        
#     botn_records.append({
#         "interactors": [
#             {'id': target_seq2id[target], 'role': 'target', 'seq': target},
#             {'id': peptide_id, 'role': 'binder', 'seq': seq},
#         ],
#         "attributes": attrs,
#         "binding": [int(x) for x in categories][0]
#     })

# print(f"[number of records]", len(botn_records))
# print(f"[number of conflict]", conflict_cases)
    
# json.dump(botn_records, open("botn.json", "w"))
    