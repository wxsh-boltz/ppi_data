import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import json
import numpy as np
from pathlib import Path

summary_file = "sabdab_summary_all.tsv"
summary_df = pd.read_csv(summary_file, sep="\t")
summary_df = summary_df.fillna("")
saving_dir = Path("processed")
saving_dir.mkdir(exist_ok=True, parents=True)
output_json_file = saving_dir / summary_file.replace(".tsv", ".json")

print(f"Read {len(summary_df)} interactions from summary.")

all_pdb_protein_seq_path = "sequences_from_raw_pdb_boltz_extended.fasta"
all_pdb_protein_pdbx_seq_path = "sequences_from_raw_pdb_pdbx_code.fasta"

# all_pdb_protein_seq_path = "sequences_from_raw_pdb_gemmi.fasta"
# all_pdb_protein_pdbx_seq_path = "sequences_from_raw_pdb_gemmi_pdbx_code.fasta"

def read_pdb_id2seq(path):
    pdb_id2seq = defaultdict(dict)
    for record in SeqIO.parse(path, "fasta"):
        pdb_id, chain_id = record.id.split("_")
        pdb_id, chain_id = pdb_id.lower(), chain_id
        if str(record.seq):
            pdb_id2seq[pdb_id.lower()][chain_id] = str(record.seq)
    
    # e.g., AA/HH etc... NOT SURE if it is correct
    for pdb_id in pdb_id2seq:
        _new_dict = {**pdb_id2seq[pdb_id]}
        for chain_id in pdb_id2seq[pdb_id]:
            seq = pdb_id2seq[pdb_id][chain_id]
            if len(chain_id) == 2: 
                for _chain_id in chain_id:
                    if _chain_id not in _new_dict and seq:
                        _new_dict[_chain_id] = seq
        pdb_id2seq[pdb_id] = _new_dict
    
    return pdb_id2seq

pdb_id2seq = read_pdb_id2seq(all_pdb_protein_seq_path)
pdb_id2pdbx_seq = read_pdb_id2seq(all_pdb_protein_pdbx_seq_path)

pdbx_seq2seq = dict()
for pdb_id in pdb_id2pdbx_seq:
    for chain_id in pdb_id2pdbx_seq[pdb_id]:
        pdbx_seq = pdb_id2pdbx_seq[pdb_id][chain_id]
        pdb_seq = pdb_id2seq[pdb_id][chain_id]
        if pdbx_seq not in pdbx_seq2seq:
            pdbx_seq2seq[pdbx_seq] = pdb_seq
        else:
            assert pdbx_seq2seq[pdbx_seq] == pdb_seq

def get_chain_seq_relax(pdb_id, chain_id, pdb_id2seq):
    if chain_id in pdb_id2seq[pdb_id]:
        return pdb_id2seq[pdb_id][chain_id]
    elif chain_id.upper() in pdb_id2seq[pdb_id]:
        return pdb_id2seq[pdb_id][chain_id.upper()]

hchain_id2seq_pdbx = dict()
hchain_pdbx_seq2id = dict()
lchain_id2seq_pdbx = dict()
lchain_pdbx_seq2id = dict()
antigen_id2seq_pdbx = dict()
antigen_pdbx_seq2id = dict()

records = defaultdict(list)
fail_chains = set()
fail_tuples = 0
for _, row in summary_df.iterrows():
    row = row.to_dict()
    
    antigen_type = row["antigen_type"]
    if "protein" not in antigen_type and "peptide" not in antigen_type:
        continue
    
    hchain = row["Hchain"]
    lchain = row["Lchain"]
    pdb = row["pdb"].lower()
    antigen_chains = row["antigen_chain"].split("|")
    antigen_types = row["antigen_type"].split("|")

    hchain_seq_pdbx, lchain_seq_pdbx = None, None

    if hchain and hchain != "NA":
        # hchain_seq = get_chain_seq_relax(pdb, hchain, pdb_id2seq) #  pdb_id2seq[pdb][hchain]
        hchain_seq_pdbx = get_chain_seq_relax(pdb, hchain, pdb_id2pdbx_seq)
        if hchain_seq_pdbx is not None:
            hchain_id2seq_pdbx[f"{pdb}_{hchain}"] = hchain_seq_pdbx # (hchain_seq, hchain_seq_pdbx)
            hchain_pdbx_seq2id[hchain_seq_pdbx] = f"{pdb}_{hchain}"
        else:
            fail_chains.add((pdb, hchain, "hchain"))
    
    if lchain and lchain != "NA":
        # lchain_seq = get_chain_seq_relax(pdb, lchain, pdb_id2seq) # pdb_id2seq[pdb][lchain]
        lchain_seq_pdbx = get_chain_seq_relax(pdb, lchain, pdb_id2pdbx_seq)
        if lchain_seq_pdbx is not None:
            lchain_id2seq_pdbx[f"{pdb}_{lchain}"] = lchain_seq_pdbx # (lchain_seq, lchain_seq_pdbx)
            lchain_pdbx_seq2id[lchain_seq_pdbx] = f"{pdb}_{lchain}"
        else:
            fail_chains.add((pdb, lchain, "lchain"))

    # antigen_seqs = []
    antigen_seqs_pdbx = []
    antigen_chains_new = []
    antigen_types_new = []

    for antigen_type, antigen_chain in zip(antigen_types, antigen_chains):
        antigen_chain = antigen_chain.strip()
        antigen_type = antigen_type.strip()
        if antigen_chain != "NA" and antigen_chain and antigen_type in ("protein", "peptide"):
            # antigen_seq = get_chain_seq_relax(pdb, antigen_chain, pdb_id2seq) # pdb_id2seq[pdb][antigen_chain]
            antigen_seq_pdbx = get_chain_seq_relax(pdb, antigen_chain, pdb_id2pdbx_seq)
            if antigen_seq_pdbx is not None:
                # antigen_seqs.append(antigen_seq)
                antigen_seqs_pdbx.append(antigen_seq_pdbx)
                antigen_chains_new.append(antigen_chain)
                antigen_types_new.append(antigen_type)

                antigen_id2seq_pdbx[f"{pdb}_{antigen_chain}"] = antigen_seq_pdbx # (antigen_seq, antigen_seq_pdbx)
                antigen_pdbx_seq2id[antigen_seq_pdbx] = f"{pdb}_{antigen_chain}"
            else:
                fail_chains.add((pdb, antigen_chain, "antigen_chain"))
    
    # row["hchain_seq"] = hchain_seq
    # row["lchain_seq"] = lchain_seq
    # row["antigen_seqs"] = antigen_seqs

    if hchain_seq_pdbx is not None and lchain_seq_pdbx is not None:
        if hchain_seq_pdbx == lchain_seq_pdbx and not row["scfv"]:
            print("Skip: lchain_seq == hchain_seq", row)
            continue
    if hchain_seq_pdbx is not None:
        if hchain_seq_pdbx in antigen_seqs_pdbx:
            print(f"Skip: hchain_seq IN antigen_seqs: {pdb}, {hchain}, {antigen_chains}")
            continue
    if lchain_seq_pdbx is not None:
        if lchain_seq_pdbx in antigen_seqs_pdbx:
            print(f"Skip: lchain_seq IN antigen_seqs: {pdb}, {lchain}, {antigen_chains}")
            continue
    
    
    sorted_indices = np.argsort(antigen_seqs_pdbx)
    antigen_seqs_pdbx = tuple(np.asarray(antigen_seqs_pdbx)[sorted_indices])
    # antigen_seqs = tuple(np.asarray(antigen_seqs)[sorted_indices])
    antigen_chains_new = np.asarray(antigen_chains_new)[sorted_indices]
    antigen_types_new = np.asarray(antigen_types_new)[sorted_indices]

    row["antigen_chain_match"] = antigen_chains_new.tolist()
    row["antigen_type_match"] = antigen_types_new.tolist()

    # if len(sorted_indices) > 1:
    #     print(sorted_indices, row["antigen_chain_match"], row["antigen_type_match"])
    #     print(row)

    if (hchain_seq_pdbx is not None or lchain_seq_pdbx is not None) and len(antigen_seqs_pdbx) >= 1:
        records[(hchain_seq_pdbx, lchain_seq_pdbx, antigen_seqs_pdbx)].append(row)
    else:
        fail_tuples += 1


print("fail_chains", len(fail_chains), fail_chains)
print(f"fail_tuples", fail_tuples)

# hchain_seq2id = dict()
# for id, seq in hchain_id2seq.items():
#     hchain_seq2id[seq] = id

with open(saving_dir / "hchain_seqs.fasta", "w") as fout, open(saving_dir / "hchain_seqs.pdbx.fasta", "w") as fout_pdbx:
    for pdbx_seq, id in hchain_pdbx_seq2id.items():
        fout_pdbx.write(f">{id}\n{pdbx_seq}\n")
        seq = pdbx_seq2seq[pdbx_seq]
        fout.write(f">{id}\n{seq}\n")

# with open("hchain_seqs.fasta", "w") as fout:
#     for seq, id in hchain_seq2id.items():
#         fout.write(f">{id}\n{seq}\n")


# lchain_seq2id = dict()
# for id, seq in lchain_id2seq.items():
#     lchain_seq2id[seq] = id
# with open("lchain_seqs.fasta", "w") as fout:
#     for seq, id in lchain_seq2id.items():
#         fout.write(f">{id}\n{seq}\n")
with open(saving_dir / "lchain_seqs.fasta", "w") as fout, open(saving_dir / "lchain_seqs.pdbx.fasta", "w") as fout_pdbx:
    for pdbx_seq, id in lchain_pdbx_seq2id.items():
        fout_pdbx.write(f">{id}\n{pdbx_seq}\n")
        seq = pdbx_seq2seq[pdbx_seq]
        fout.write(f">{id}\n{seq}\n")

# antigen_seq2id = dict()
# for id, seq in antigen_id2seq.items():
#     antigen_seq2id[seq] = id
# with open("antigen_seqs.fasta", "w") as fout:
#     for seq, id in antigen_seq2id.items():
#         fout.write(f">{id}\n{seq}\n")
with open(saving_dir / "antigen_seqs.fasta", "w") as fout, open(saving_dir / "antigen_seqs.pdbx.fasta", "w") as fout_pdbx:
    for pdbx_seq, id in antigen_pdbx_seq2id.items():
        fout_pdbx.write(f">{id}\n{pdbx_seq}\n")
        seq = pdbx_seq2seq[pdbx_seq]
        fout.write(f">{id}\n{seq}\n")

print(f"Interactions with protein/peptide antigens: {len(records)}")
print(f"Hchains (unique ID): {len(hchain_id2seq_pdbx)}, Hchains (unique seqs): {len(hchain_pdbx_seq2id)}")
print(f"Lchains (unique ID): {len(lchain_id2seq_pdbx)}, Lchains (unique seqs): {len(lchain_pdbx_seq2id)}")
print(f"Antigen (unique ID): {len(antigen_id2seq_pdbx)}, Antigen (unique seqs): {len(antigen_pdbx_seq2id)}")

output_records = []
for (hchain_pdbx, lchain_pdbx, antigens_pdbx), attributes in records.items():
    hchain = pdbx_seq2seq[hchain_pdbx] if hchain_pdbx is not None else None
    lchain = pdbx_seq2seq[lchain_pdbx] if lchain_pdbx is not None else None
    # hchain, lchain = pdbx_seq2seq[hchain_pdbx], pdbx_seq2seq[lchain_pdbx]

    if hchain_pdbx != lchain_pdbx:
        interactors = []
        if hchain_pdbx is not None:
            interactors.append({"seq": hchain, "id": hchain_pdbx_seq2id[hchain_pdbx], "role": "hchain", "seq_pdbx": hchain_pdbx})
        if lchain is not None:
            interactors.append({"seq": lchain, "id": lchain_pdbx_seq2id[lchain_pdbx], "role": "lchain", "seq_pdbx": lchain_pdbx})
        for antigen_pdbx in antigens_pdbx:
            antigen = pdbx_seq2seq[antigen_pdbx]
            interactors.append({"seq": antigen, "id": antigen_pdbx_seq2id[antigen_pdbx], "role": "antigen", "seq_pdbx": antigen_pdbx})
        
        output_records.append({"interactors": interactors, "attributes": attributes})
    else: # scfv: hchian and lchain are in the same chain, connected by a flexible linker
        interactors = []
        for attr in attributes:
            assert attr["scfv"] # double check
        interactors.append({"seq": hchain, "id": hchain_pdbx_seq2id[hchain_pdbx], "role": "fv", "seq_pdbx": hchain_pdbx})
        for antigen_pdbx in antigens_pdbx:
            antigen = pdbx_seq2seq[antigen_pdbx]
            interactors.append({"seq": antigen, "id": antigen_pdbx_seq2id[antigen_pdbx], "role": "antigen", "seq_pdbx": antigen_pdbx})
        output_records.append({"interactors": interactors, "attributes": attributes})

print(len(output_records))
print("H/L pairs:", len(set([(x[0], x[1]) for x in records])))
json.dump(output_records, open(output_json_file, "w"))