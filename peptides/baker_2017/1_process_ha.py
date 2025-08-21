import pandas as pd
from collections import Counter, defaultdict
from Bio import SeqIO
import json
from pathlib import Path
from Bio.PDB import PDBParser, PPBuilder
from tqdm import tqdm
import gemmi

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

def extract_protein_sequences(pdb_file):
    """
    Extract protein sequences for each chain from a PDB file using Gemmi.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        dict: Dictionary with chain IDs as keys and sequences as values
    """
    # Read the PDB file
    if isinstance(pdb_file, Path):
        pdb_file = str(pdb_file)

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



debug = False

root_dir = Path("/data/scratch/wxsh/data/peptides/nature_2017_baker/archive_MassivelyParallelDeNovoProteinDesignForTargetedTherapeutics_2017/")
path_ha_gen1 = root_dir / "03_HighThroughputDesign_FluH1_Gen1_relax.dat"
path_ha_gen2 = root_dir / "06_HighThroughputDesign_FluH1_Gen2_relax.dat"
csv_path = Path("ha.csv")
dataset_name = "ha"

if not csv_path.exists():
    df_ha_gen1 = parse(path_ha_gen1)
    df_ha_gen2 = parse(path_ha_gen2)
    print("HA_GEN1", len(df_ha_gen1), "HA_GEN2", len(df_ha_gen2))
    print("HA_GEN1 seq set:", len(set(df_ha_gen1["Sequence"])), "HA_GEN2 seq set:", len(set(df_ha_gen2["Sequence"])))
    print(Counter(df_ha_gen1["Category"]).most_common(), Counter(df_ha_gen2["Category"]).most_common())
    if debug:
        df_ha_gen1 = df_ha_gen1[:100]
        df_ha_gen2 = df_ha_gen2[:100]

    ha_gen1_pdb_dir = root_dir / Path("04_HighThroughputDesign_FluH1_Gen1_relax_models")
    ha_gen2_pdb_dir = root_dir / Path("07_HighThroughputDesign_FluH1_Gen2_relax_models")
    ha_gen1_pdb_files = list(ha_gen1_pdb_dir.rglob("*.pdb"))
    ha_gen2_pdb_files = list(ha_gen2_pdb_dir.rglob("*.pdb"))
    print("HA GEN1 pdb_files:", len(ha_gen1_pdb_files), "HA GEN2 pdb_files:", len(ha_gen2_pdb_files))

    ha_id2seq_pairs = dict()
    for id in tqdm(df_ha_gen1["#ID"]):
        pdb_file = ha_gen1_pdb_dir / (id + ".pdb")
        pdb_file_nonrelaxed = ha_gen1_pdb_dir / (id + "_nonrelaxed.pdb") # HB1.4020.0_nonrelaxed.pdb
        if not pdb_file.exists() and not pdb_file_nonrelaxed.exists():
            print("File not exists:", pdb_file)
            continue
        if pdb_file.exists():
            sequences = extract_protein_sequences(pdb_file)
        else:
            sequences = extract_protein_sequences(pdb_file_nonrelaxed)
        assert len(sequences) == 2
        if id in ha_id2seq_pairs:
            assert ha_id2seq_pairs[id] == (sequences["A"], sequences["B"])
        ha_id2seq_pairs[id] = (sequences["A"], sequences["B"])

    for id in tqdm(df_ha_gen2["#ID"]):
        pdb_file = ha_gen2_pdb_dir / (id + ".pdb")
        if not pdb_file.exists():
            print("File not exists:", pdb_file)
            continue
        sequences = extract_protein_sequences(pdb_file)
        assert len(sequences) == 2
        if id in ha_id2seq_pairs:
            assert ha_id2seq_pairs[id] == (sequences["A"], sequences["B"])
        ha_id2seq_pairs[id] = (sequences["A"], sequences["B"])

    print(f"Read {len(ha_id2seq_pairs)} seq pairs from PDB")

    # for i, pdb_file in enumerate(tqdm(ha_gen1_pdb_files + ha_gen2_pdb_files)):
    #     if debug and i == 1000: # debug
    #         break

    #     id = pdb_file.stem
    #     sequences = extract_protein_sequences(pdb_file)
    #     assert len(sequences) == 2
    #     ha_id2seq_pairs[id] = (sequences["A"], sequences["B"])

    df_ha = pd.concat([df_ha_gen1, df_ha_gen2], axis=0)
    print("Merged G1+G2:", len(df_ha))

    target_seqs = []
    sequence_from_pdb = []
    for seq, cat, id in zip(df_ha["Sequence"], df_ha["Category"], df_ha["#ID"]):
        if id not in ha_id2seq_pairs:
            print(id)
            continue
        target, peptide = ha_id2seq_pairs[id]
        target_seqs.append(target)
        sequence_from_pdb.append(peptide)
    df_ha["target_seq"] = target_seqs
    df_ha["sequence_from_pdb"] = sequence_from_pdb
    df_ha.to_csv(csv_path, index=False)
else:
    df_ha = pd.read_csv(csv_path)

print(f"len(df_ha): {len(df_ha)}")

# print(f"[Unique seqs] num of binders: {len(set(df_ha['Sequence']))}, num of targets {len(set(df_ha['target_seq']))}, num of binders from pdb {len(set(df_ha['sequence_from_pdb']))}")


# seq2categories = defaultdict(list)
# peptide_seq2id = dict()
# all_target_seqs = set()
# for seq, cat, id, target, sequence_from_pdb in zip(df_ha["Sequence"], df_ha["Category"], df_ha["#ID"], df_ha["target_seq"], df_ha["sequence_from_pdb"]):
#     if seq not in peptide_seq2id:
#         peptide_seq2id[seq] = id
#     seq2categories[(seq, target)].append((cat, id, sequence_from_pdb))
#     all_target_seqs.add(target)
# print(f"[Number of seq pairs]: {len(seq2categories)}")

# count = dict(Counter(df_ha["Category"])) # .most_common()
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


# records = []
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
        
#     records.append({
#         "interactors": [
#             {'id': target_seq2id[target], 'role': 'target', 'seq': target},
#             {'id': peptide_id, 'role': 'binder', 'seq': seq},
#         ],
#         "attributes": attrs,
#         "binding": [int(x) for x in categories][0]
#     })

# print(f"[number of records]", len(records))
# print(f"[number of conflict]", conflict_cases)
    
# json.dump(records, open(f"{dataset_name}.json", "w"))



# seq2categories = defaultdict(list)
# seq2id = dict()
# for seq, cat, id, target, sequence_from_pdb in zip(df_ha["Sequence"], df_ha["Category"], df_ha["#ID"], df_ha["target_seq"], df_ha["sequence_from_pdb"]):
#     if seq not in seq2id:
#         seq2id[seq] = id
#     seq2categories[(seq, target)].append((cat, id, sequence_from_pdb))

# print("Seq pairs:", len(seq2categories))

# records = []
# conflict_cases = 0
# for seq, target in seq2categories:
#     peptide_id = seq2id[seq]
    
#     attrs = []
#     categories = []
#     for cat, id, sequence_from_pdb in seq2categories[(seq, target)]:
#         categories.append(cat)
#         attrs.append({"category": cat, "id": id, "sequence_from_pdb": sequence_from_pdb})

#     categories_set = set(categories)
#     if len(categories_set) > 1:
#         print(categories_set)
#         conflict_cases += 1
    
#     records.append({
#         "interactors": [
#             {'id': 'botn', 'role': 'target', 'seq': target},
#             {'id': peptide_id, 'role': 'binder', 'seq': seq},
#         ],
#         "attributes": attrs,
#         "binding": [int(x) for x in categories]
#     })

# print(len(records), conflict_cases)

# json.dump(records, open("ha.json", "w"))