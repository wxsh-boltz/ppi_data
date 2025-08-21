from collections import defaultdict
import json, argparse
from tqdm import tqdm
import numpy as np
from Bio import SeqIO
from pathlib import Path

def read_alignment(path, identity_threshold=0.7, skip_self=False):
    query_id2similar_target_ids = defaultdict(set)
    evalues = []
    with open(path) as fin:
        for line in fin:
            query_id, target_id, identity = line.strip().split()[:3]
            
            identity = float(identity)
            if identity < identity_threshold:
                continue
            if skip_self and query_id == target_id:
                continue
            
            if len(line.strip().split()) > 10:
                evalue = float(line.strip().split()[10])
                if evalue > 0.1:
                    continue
                evalues.append(evalue)

            query_id2similar_target_ids[query_id].add(target_id)
    if len(evalues) > 0:
        print(f"Evalue: max {np.max(evalues)}, mean {np.mean(evalues)}")
    return query_id2similar_target_ids

# def quality_control_seq(seq, min_seq_len=30, max_X_ratio = 0.5, replace_non_standard_aa_by_X = True, exclude_set=None):
#     if len(seq) < min_seq_len:
#         return False
#     if exclude_set is not None and seq in exclude_set:
#         return False
#     if replace_non_standard_aa_by_X:
#         seq_new = []
#         for c in seq:
#             if c.islower():
#                 seq_new.append("X")
#             else:
#                 seq_new.append(c)
#         seq = "".join(seq_new)
        
#     if sum(np.asarray(list(seq)) == "X") >= max_X_ratio:
#         return False
#     # if sum(np.char.islower(list(seq))) >= max_X_ratio:
#         # return False

#     return True 

def quality_control_seq(seq, max_X_ratio=0.5, replace_non_standard_aa_by_X=False):
    if replace_non_standard_aa_by_X:
        seq_new = []
        for c in seq:
            if c.islower():
                seq_new.append("X")
            else:
                seq_new.append(c)
        seq = "".join(seq_new)
    
    if "*" in seq.rstrip("*"):
        return False

        
    if np.mean(np.asarray(list(seq)) == "X") >= max_X_ratio:
        return False
    # if sum(np.char.islower(list(seq))) >= max_X_ratio:
        # return False

    return True 

def entry_id2chain_ids(rcsb_ids):
    # rcsb_ids: List of string, like xxxx_A, xxxx_B
    r = defaultdict(list)
    for id in rcsb_ids:
        r[id.split("_")[0]].append(id.split("_")[1])
    return r


parser = argparse.ArgumentParser()
parser.add_argument('--input_path', type=str)
parser.add_argument('--rcsb_ids_path', type=str, default="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/rcsb_boltz2_train_ids.txt")
parser.add_argument('--output_path', type=str, default=None)
parser.add_argument('--vh_aln_paths', type=str, nargs='+')
parser.add_argument('--vl_aln_paths', type=str, default=None, nargs='+')
parser.add_argument('--antigen_aln_paths', type=str, nargs='+')
parser.add_argument('--vh_self_aln_path', type=str, default=None)
parser.add_argument('--vl_self_aln_path', type=str, default=None)
parser.add_argument('--antigen_self_aln_path', type=str, default=None)
parser.add_argument('--vh_sequence_identity_threshold', type=float, default=0.7)
parser.add_argument('--vl_sequence_identity_threshold', type=float, default=0.7)
parser.add_argument('--antigen_sequence_identity_threshold', type=float, default=0.7)
### TODO: special process for the shorter antigens?
parser.add_argument('--peptide_antigen_aln_path', type=str, default=None)
parser.add_argument('--peptide_antigen_list_path', type=str, default=None)
parser.add_argument('--assay_id', type=str)

args = parser.parse_args()

rcsb_ids = set()
with open(args.rcsb_ids_path) as fin:
    for line in fin:
        rcsb_ids.add(line.strip().lower())
print(f"RCSB ids: {len(rcsb_ids)}")

antigen_id2similar_rcsb_ids = defaultdict(set)
for antigen_aln_path in args.antigen_aln_paths:
    _res = read_alignment(antigen_aln_path, identity_threshold=args.antigen_sequence_identity_threshold)
    for key in _res:
        antigen_id2similar_rcsb_ids[key].update(_res[key])
print(f"Read {len(antigen_id2similar_rcsb_ids)} antigen alignments")

if args.peptide_antigen_aln_path is not None and args.peptide_antigen_list_path:
    peptide_antigen_ids = set()
    with open(args.peptide_antigen_list_path) as fin:
        for line in fin:
            id = line.strip().split()[0]
            peptide_antigen_ids.add(id)
    peptide_antigen_id2similar_rcsb_ids = read_alignment(args.peptide_antigen_aln_path, identity_threshold=-100)
    print(f"Read {len(peptide_antigen_id2similar_rcsb_ids)} peptide antigen alignments")
else:
    peptide_antigen_id2similar_rcsb_ids, peptide_antigen_ids = None, None

vh_id2similar_rcsb_ids = defaultdict(set)
for vh_aln_path in args.vh_aln_paths:
    # assert Path(vh_aln_path).exists()
    _res = read_alignment(vh_aln_path, identity_threshold=args.vh_sequence_identity_threshold)
    for key in _res:
        vh_id2similar_rcsb_ids[key].update(_res[key])
print(f"Read {len(vh_id2similar_rcsb_ids)} VH alignments")
if args.vl_aln_paths is not None: # for VHH
    vl_id2similar_rcsb_ids = defaultdict(set)
    for vl_aln_path in args.vl_aln_paths:
        _res = read_alignment(vl_aln_path, identity_threshold=args.vl_sequence_identity_threshold)
        for key in _res:
            vl_id2similar_rcsb_ids[key].update(_res[key])
    # vl_id2similar_rcsb_ids = read_alignment(args.vl_aln_path, identity_threshold=args.vl_sequence_identity_threshold)
    print(f"Read {len(vl_id2similar_rcsb_ids)} VL alignments")

if args.antigen_self_aln_path is not None:
    antigen_id2antigen_ids = read_alignment(args.antigen_self_aln_path, identity_threshold=1.0)
    print(f"Read {len(antigen_id2antigen_ids)} antigen alignments (to self)")
else:
    antigen_id2antigen_ids = None
if args.vh_self_aln_path is not None:
    vh_id2vh_ids = read_alignment(args.vh_self_aln_path, identity_threshold=1.0)
    print(f"Read {len(vh_id2vh_ids)} VH alignments (to self)")
else:
    vh_id2vh_ids = None
if args.vl_self_aln_path is not None:
    vl_id2vl_ids = read_alignment(args.vl_self_aln_path, identity_threshold=1.0)
    print(f"Read {len(vl_id2vl_ids)} VL alignments (to self)")
else:
    vl_id2vl_ids = None


input_data = json.load(open(args.input_path))
print(f"Read {len(input_data)} input data")
saving_dir = Path(args.input_path.replace(".json", f".test_vh_{args.vh_sequence_identity_threshold}_vl_{args.vl_sequence_identity_threshold}_an_{args.antigen_sequence_identity_threshold}"))
saving_dir.mkdir(exist_ok=True, parents=True)
if args.output_path is None:
    saving_path = saving_dir / "data.json"
else:
    saving_path = args.output_path


data_not_in_pdb = []
antigen_skip_by_mmseqs = set()
vh_skip_by_mmseqs = set()
vl_skip_by_mmseqs = set()
antigen_skip_by_qc = set()
vh_skip_by_qc = set()
vl_skip_by_qc = set()
pos_data = []
neg_data = []
for data in tqdm(input_data):
    ### 1. exclude by pdb ids:
    skip_pdb = False
    if "attributes" in data:
        for attribute in data["attributes"]:
            if "pdb" in attribute: # sabdab
                pdb = attribute["pdb"].lower()
                if pdb in rcsb_ids:
                    skip_pdb = True
                    break

        # elif "PDB_ID" in attribute: # snac
        #     pdb = attribute["PDB_ID"].lower()
        #     if pdb in rcsb_ids:
        #         skip_pdb = True
        #         break
    if skip_pdb:
        continue
    
    vh_ids = []
    vh_seqs = []
    vl_ids = []
    vl_seqs = []
    antigen_ids = []
    antigen_seqs = []
    for interactor in data["interactors"]:
        if interactor["role"] == "vh":
            vh_ids.append(interactor["id"])
            vh_seqs.append(interactor["seq"])
        elif interactor["role"] == "vl":
            vl_ids.append(interactor["id"])
            vl_seqs.append(interactor["seq"])
        elif interactor["role"] == "antigen":
            antigen_ids.append(interactor["id"])
            antigen_seqs.append(interactor["seq"])
    
    if "binding" in data:
        if data["binding"]:
            pos_data.append(data)
        else:
            neg_data.append(data)

    assert len(vh_ids) <= 1, vh_ids
    assert len(vl_ids) <= 1, vl_ids
    assert len(vh_ids) == 1 or len(vl_ids) == 1, data
    assert len(antigen_ids) >= 1

    if antigen_id2antigen_ids:
        _antigen_ids_skip_by_mmseqs = [(id, seq) for id, seq in zip(antigen_ids, antigen_seqs) if len(antigen_id2antigen_ids[id]) == 0]
    else:
        _antigen_ids_skip_by_mmseqs = []
    if peptide_antigen_ids is not None:
        _antigen_ids_skip_by_mmseqs = [(id, seq) for id, seq in _antigen_ids_skip_by_mmseqs if id not in peptide_antigen_ids]

    if vh_id2vh_ids:
        _vh_ids_skip_by_mmseqs = [(id, seq) for id, seq in zip(vh_ids, vh_seqs) if len(vh_id2vh_ids[id]) == 0]
    else:
        _vh_ids_skip_by_mmseqs = []
    if vl_id2vl_ids:
        _vl_ids_skip_by_mmseqs = [(id, seq) for id, seq in zip(vl_ids, vl_seqs) if len(vl_id2vl_ids[id]) == 0]
    else:
        _vl_ids_skip_by_mmseqs = []

    skip_mmseqs = False
    if len(_antigen_ids_skip_by_mmseqs) > 0:
        antigen_skip_by_mmseqs.update(_antigen_ids_skip_by_mmseqs)
        skip_mmseqs = True
    if len(_vh_ids_skip_by_mmseqs) > 0:
        vh_skip_by_mmseqs.update(_vh_ids_skip_by_mmseqs)
        skip_mmseqs = True
    if len(_vl_ids_skip_by_mmseqs) > 0:
        vl_skip_by_mmseqs.update(_vl_ids_skip_by_mmseqs)
        skip_mmseqs = True
    
    if skip_mmseqs:
        continue

    if any(not quality_control_seq(seq) for seq in antigen_seqs) or any(not quality_control_seq(seq) for seq in vh_seqs) or any(not quality_control_seq(seq) for seq in vl_seqs):
        antigen_skip_by_qc.update(seq for seq in antigen_seqs if not quality_control_seq(seq))
        vh_skip_by_qc.update(seq for seq in vh_seqs if not quality_control_seq(seq))
        vl_skip_by_qc.update(seq for seq in vl_seqs if not quality_control_seq(seq))
        continue

    if len(vh_ids) > 0:
        vh_rcsb_ids = vh_id2similar_rcsb_ids.get(vh_ids[0], set())
    else:
        vh_rcsb_ids = set()
    vh_rcsb_ids = entry_id2chain_ids(vh_rcsb_ids)
    if len(vl_ids) > 0:
        vl_rcsb_ids = vl_id2similar_rcsb_ids.get(vl_ids[0], set())
    else:
        vl_rcsb_ids = set()
    vl_rcsb_ids = entry_id2chain_ids(vl_rcsb_ids)
    
    antigen_rcsb_ids = set()
    for antigen_id in antigen_ids:
        antigen_rcsb_ids.update(antigen_id2similar_rcsb_ids.get(antigen_id, set()))
        if peptide_antigen_id2similar_rcsb_ids is not None:
            # if len(antigen_id2similar_rcsb_ids.get(antigen_id, set())) > 0:
            #     if len(peptide_antigen_id2similar_rcsb_ids.get(antigen_id, set())) > 0:
            #         print(antigen_id)
            #         print(peptide_antigen_id2similar_rcsb_ids.get(antigen_id, set()) - antigen_id2similar_rcsb_ids.get(antigen_id, set()))
            #         print(set(vh_rcsb_ids.keys()))
            #         print(set(vl_rcsb_ids.keys()))
            #         exit()
            antigen_rcsb_ids.update(peptide_antigen_id2similar_rcsb_ids.get(antigen_id, set()))
    antigen_rcsb_ids = entry_id2chain_ids(antigen_rcsb_ids)
        
    if len(vl_ids) != 0 and len(vh_ids) != 0: # for antibody, both VH and VL exist
        overlap_entries = set(vh_rcsb_ids.keys()) & set(vl_rcsb_ids.keys()) & set(antigen_rcsb_ids.keys())
    elif len(vl_ids) == 0: # nanobody with VH only
        overlap_entries = set(vh_rcsb_ids.keys()) & set(antigen_rcsb_ids.keys())
    elif len(vh_ids) == 0: # engineered antibody with VL only
        overlap_entries = set(vl_rcsb_ids.keys()) & set(antigen_rcsb_ids.keys())
    
    if len(overlap_entries) == 0:
        data["assay_id"] = args.assay_id
        data_not_in_pdb.append(data)

print(f"pos: {len(pos_data)}, neg: {len(neg_data)}")

def print_skip_by_mmseqs_to_file(path, skip_by_mmseqs):
    if len(skip_by_mmseqs) > 0:
        with open(path, "w") as fout:
            for id, s in skip_by_mmseqs:
                fout.write(f"{id}\t{s}\n")

print(f"antigen_skip_by_mmseqs: {len(antigen_skip_by_mmseqs)}")
print_skip_by_mmseqs_to_file(Path(args.antigen_aln_paths[0]).parent / "antigen_skip_by_mmseqs.txt", antigen_skip_by_mmseqs)
print(f"vh_skip_by_mmseqs: {len(vh_skip_by_mmseqs)}")
print_skip_by_mmseqs_to_file(Path(args.vh_aln_paths[0]).parent / "vh_skip_by_mmseqs.txt", vh_skip_by_mmseqs)
print(f"vl_skip_by_mmseqs: {len(vl_skip_by_mmseqs)}")
if args.vl_aln_paths is not None and len(args.vl_aln_paths) > 0:
    print_skip_by_mmseqs_to_file(Path(args.vl_aln_paths[0]).parent / "vl_skip_by_mmseqs.txt", vl_skip_by_mmseqs)
print(f"antigen_skip_by_qc: {len(antigen_skip_by_qc)}")
print(f"vh_skip_by_qc: {len(vh_skip_by_qc)}")
print(f"vl_skip_by_qc: {len(vl_skip_by_qc)}")

print(f"data_not_in_pdb: {len(data_not_in_pdb)}")
if len(data_not_in_pdb) > 0 and "binding" in data_not_in_pdb[0]:
    new_pos_num = len([x for x in data_not_in_pdb if x['binding']])
    new_neg_num = len([x for x in data_not_in_pdb if not x['binding']])
    print(f"pos: {new_pos_num} ({new_pos_num/len(pos_data)}), neg: {new_neg_num} ({new_neg_num/len(neg_data)})")

json.dump(data_not_in_pdb, open(saving_path, "w"))

test_antigens = dict()
test_vhs = dict()
test_vls = dict()
for data in data_not_in_pdb:
    for interactor in data["interactors"]:
        if interactor["role"] == "vh":
            test_vhs[interactor["id"]] = interactor["seq"]
        elif interactor["role"] == "vl":
            test_vls[interactor["id"]] = interactor["seq"]
        elif interactor["role"] == "antigen":
            test_antigens[interactor["id"]] = interactor["seq"]
print(f"Test antigens: {len(test_antigens)}")
print(f"Test vhs: {len(test_vhs)}")
print(f"Test vls: {len(test_vls)}")
with open(saving_dir / "vh_seqs.fasta", "w") as fout:
    for id, seq in test_vhs.items():
        fout.write(f">{id}\n{seq}\n")
with open(saving_dir / "vl_seqs.fasta", "w") as fout:
    for id, seq in test_vls.items():
        fout.write(f">{id}\n{seq}\n")
with open(saving_dir / "antigen_seqs.fasta", "w") as fout:
    for id, seq in test_antigens.items():
        fout.write(f">{id}\n{seq}\n")



## subset of novel antigens?
rcsb_training_multimer_id_path = "/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/rcsb_boltz2_train_multimers_id.txt"
rcsb_training_multimer_ids = set()
with open(rcsb_training_multimer_id_path) as fin:
    for line in fin:
        rcsb_training_multimer_ids.update(line.strip().split())

novel_antigens = set()
for antigen_id in test_antigens:
    antigen_rcsb_ids = antigen_id2similar_rcsb_ids.get(antigen_id, set())
    if peptide_antigen_id2similar_rcsb_ids is not None:
        antigen_rcsb_ids.update(peptide_antigen_id2similar_rcsb_ids.get(antigen_id, set()))
    if len(set(antigen_rcsb_ids) & rcsb_training_multimer_ids) == 0:
        novel_antigens.add(antigen_id)
        
print(f"Novel antigens? {len(novel_antigens)}")
with open(saving_dir / "antigen_seqs_novel.fasta", "w") as fout:
    for id, seq in test_antigens.items():
        if id in novel_antigens:
            fout.write(f">{id}\n{seq}\n")
