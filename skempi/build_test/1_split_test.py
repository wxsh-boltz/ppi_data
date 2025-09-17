from collections import defaultdict
import json, argparse
from tqdm import tqdm
import numpy as np
from Bio import SeqIO
from pathlib import Path

def read_alignment(path, identity_threshold=0.7, skip_self=False):
    query_id2similar_target_ids = defaultdict(set)
    with open(path) as fin:
        for line in fin:
            query_id, target_id, identity = line.strip().split()[:3]
            identity = float(identity)
            if identity < identity_threshold:
                continue
            if skip_self and query_id == target_id:
                continue
            query_id2similar_target_ids[query_id].add(target_id)
    return query_id2similar_target_ids

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

    return True 

def entry_id2chain_ids(rcsb_ids):
    # rcsb_ids: List of string, like xxxx_A, xxxx_B
    r = defaultdict(list)
    for id in rcsb_ids:
        r[id.split("_")[0]].append(id.split("_")[1])
    return r

parser = argparse.ArgumentParser()
parser.add_argument('--input_path', type=str)
parser.add_argument('--output_path', type=str, default=None)
parser.add_argument('--aln_paths', type=str, nargs='+')
parser.add_argument('--self_aln_path', type=str, default=None)
parser.add_argument('--peptide_aln_path', type=str, default=None)
parser.add_argument('--peptide_path', type=str, default=None)
parser.add_argument('--sequence_identity_threshold', type=float, default=0.7)

args = parser.parse_args()

protein_id2similar_rcsb_ids = defaultdict(set)
for aln_path in args.aln_paths:
    _res = read_alignment(aln_path, identity_threshold=args.sequence_identity_threshold)
    for key in _res:
        protein_id2similar_rcsb_ids[key].update(_res[key])
print(f"Read {len(protein_id2similar_rcsb_ids)} alignments")
if args.peptide_aln_path is not None:
    pep_id2similar_rcsb_ids = read_alignment(args.peptide_aln_path, identity_threshold=-100)
    print(f"Read {len(pep_id2similar_rcsb_ids)} peptide alignments")
else:
    pep_id2similar_rcsb_ids = None
if args.peptide_path is not None:
    peptides = set([x.strip() for x in open(args.peptide_path)])
else:
    peptides = None
    
if args.self_aln_path is not None:
    protein_id2self_ids = read_alignment(args.self_aln_path, identity_threshold=1.0)
    print(f"Read {len(protein_id2self_ids)} alignments (to self)")
else:
    protein_id2self_ids = None

input_data = json.load(open(args.input_path))
print(f"Read {len(input_data)} input data")
if args.output_path is None:
    saving_dir = Path(args.input_path.replace(".json", f".test_si_{args.sequence_identity_threshold}"))
    saving_dir.mkdir(exist_ok=True, parents=True)
    # saving_path = args.input_path.replace(".json", f".test_si_{args.sequence_identity_threshold}.pos.json")
    saving_path = saving_dir / "data.json"
else:
    saving_dir = args.output_path.parent
    saving_path = args.output_path



data_not_in_pdb = []
skip_by_mmseqs = set()
skip_by_qc = set()
pos_data = []
neg_data = []
for data in tqdm(input_data):
    _ids = []
    _seqs = []
    for interactor in data["interactors"]:
        _ids.append(interactor["id"])
        _seqs.append(interactor["seq"])

    if protein_id2self_ids:
        _ids_skip_by_mmseqs = [(id, seq) for id, seq in zip(_ids, _seqs) if id not in protein_id2self_ids]
    else:
        _ids_skip_by_mmseqs = []
    if peptides is not None:
        _ids_skip_by_mmseqs = [(id, seq) for id, seq in _ids_skip_by_mmseqs if id not in peptides]

    skip_mmseqs = False
    if len(_ids_skip_by_mmseqs) > 0:
        skip_by_mmseqs.update(_ids_skip_by_mmseqs)
        continue

    if any(not quality_control_seq(seq) for seq in _seqs):
        skip_by_qc.update(seq for seq in _seqs if not quality_control_seq(seq))
        continue

    interactor_similar_rcsb_ids = []
    for id in _ids:
        _rcsb_ids = protein_id2similar_rcsb_ids.get(id, set())
        if pep_id2similar_rcsb_ids is not None:
            _rcsb_ids.update(pep_id2similar_rcsb_ids.get(id, set()))
        interactor_similar_rcsb_ids.append(set(x.split("_")[0] for x in _rcsb_ids))

    # _not_in_pdb = True
    # for i in range(len(interactor_similar_rcsb_ids)):
    #     for j in range(i+1, len(interactor_similar_rcsb_ids)):
    #         if len(interactor_similar_rcsb_ids[i] & interactor_similar_rcsb_ids[j]) > 0:
    #             _not_in_pdb = False

    overlap_entries = set.intersection(*interactor_similar_rcsb_ids)
    if len(overlap_entries) == 0:
        data_not_in_pdb.append(data)
    # if _not_in_pdb:
        # data_not_in_pdb.append(data)

print(f"skip_by_mmseqs: {len(skip_by_mmseqs)}")
print(f"skip_by_qc: {len(skip_by_qc)}")

with open(saving_dir / "skip_by_mmseqs.log", "w") as fout:
    for id, seq in skip_by_mmseqs:
        fout.write(f"{id}\t{seq}\n")
with open(saving_dir / "skip_by_qc.log", "w") as fout:
    for seq in skip_by_qc:
        fout.write(f"{seq}\n")

print(f"data_not_in_pdb: {len(data_not_in_pdb)}")
# if len(data_not_in_pdb) > 0 and "binding" in data_not_in_pdb[0]:
#     new_pos_num = len([x for x in data_not_in_pdb if x['binding']])
#     new_neg_num = len([x for x in data_not_in_pdb if not x['binding']])
#     print(f"pos: {new_pos_num} ({new_pos_num/len(pos_data)}), neg: {new_neg_num} ({new_neg_num/len(neg_data)})")

json.dump(data_not_in_pdb, open(saving_path, "w"))

test_id2seq = dict()
for data in data_not_in_pdb:
    interactors = data["interactors"]
    for interactor in interactors:
        test_id2seq[interactor["id"]] = interactor["seq"]
with open(saving_dir / "seqs.fasta", "w") as fout:
    for id, seq in test_id2seq.items():
        fout.write(f">{id}\n{seq}\n")