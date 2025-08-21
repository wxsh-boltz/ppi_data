from typing import List
from hydra.utils import instantiate
from omegaconf import OmegaConf
import hydra, json
from typing import Any
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

def read_rcsb_complexes(rcsb_multimer_ids_path="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/rcsb_boltz2_train_multimers_id.txt"):
    complex_chain_ids = set()
    with open(rcsb_multimer_ids_path) as fin:
        for line in fin:
            chain_ids = line.strip().split()
            pdb_id = chain_ids[0].split("_")[0]
            complex_chain_ids.update(chain_ids)
            assert all(chain_id.split("_")[0] == pdb_id for chain_id in chain_ids)
            # pdb_id2chain_ids[pdb_id] = chain_ids
    return complex_chain_ids

def read_alignment(paths, identity_threshold=0.5, skip_self=False):
    query_id2similar_target_ids = defaultdict(list)
    if isinstance(paths, str) or isinstance(paths, Path):
        paths = [paths]
    for path in paths:
        with open(path) as fin:
            for line in fin:
                query_id, target_id, identity = line.strip().split()[:3]
                identity = float(identity)
                if identity < identity_threshold:
                    continue
                if skip_self and query_id == target_id:
                    continue
                query_id2similar_target_ids[query_id].append(target_id)
    return query_id2similar_target_ids

def read_target_seqs(path):
    id2seq = dict()
    for record in SeqIO.parse(path, "fasta"):
        id2seq[record.id] = str(record.seq)
    return id2seq

def get_target_ids(dataset, role):
    target_id2seq = dict()
    for data in dataset:
        for interactor in data["interactors"]:
            if interactor["role"] == role:
                if interactor["id"] in target_id2seq:
                    assert target_id2seq[interactor["id"]] == interactor["seq"]
                target_id2seq[interactor["id"]] = interactor["seq"]
    return target_id2seq

def main():
    rcsb_train_complex_chain_ids = read_rcsb_complexes()
    print(f"RCSB complex ids: {len(rcsb_train_complex_chain_ids)}, {list(rcsb_train_complex_chain_ids)[:3]}")
    cfg = OmegaConf.load("datasets.yaml")  # Just load config.yaml
    # Instantiate each dataset
    test_datasets = [instantiate(d) for d in cfg.test_datasets]
    print(f"Read {len(test_datasets)} test_datasets")
    new_target_saving_dir = Path("/data/rbg/users/wxsh/esm3/data/finalize/new_targets")
    
    for test_dataset in test_datasets:
        if test_dataset.dataset_type not in ("miniprotein", "antibody", "tcr_phmc"):
            continue
        ### To discovery new targets:
        #### (1) Normal sequences: not exclude by mmseqs, and have more than 0.7 similarity with the existing complexes;
        #### (2) Short sequences (peptides): in peptide list, have less than 3 edit distance with existing complexes.

        data = json.load(open(test_dataset.json_path))
        print(f"Process {test_dataset.dataset_name}: read {len(data)} data")
        target_id2seq = get_target_ids(data, test_dataset.target_role)
        print(f"Number of targets: {len(target_id2seq)}")

        target_id2similar_rcsb_ids = read_alignment(test_dataset.alignment_paths.target_to_rcsb)
        target_id2target_ids = read_alignment(test_dataset.alignment_paths.target_to_target)

        if test_dataset.target_peptide_seq_path is not None:
            target_peptide_ids = set([line.strip() for line in open(test_dataset.target_peptide_seq_path)])
            target_peptide_id2similar_rcsb_ids = read_alignment(test_dataset.alignment_paths.target_peptide_to_rcsb, -1)
            print("target_peptide", len(target_peptide_id2similar_rcsb_ids))
        print(len(target_id2similar_rcsb_ids), len(target_id2target_ids))
        
        novel_target_ids = set()
        for target in target_id2seq:
            assert target in target_id2target_ids or (target in target_peptide_ids if test_dataset.target_peptide_seq_path is not None else False)
            
            set_similar_rcsb_ids = set()
            if target in target_id2target_ids: # not exclude by mmseqs
                set_similar_rcsb_ids.update(target_id2similar_rcsb_ids[target])
            
            if test_dataset.target_peptide_seq_path is not None and target in target_peptide_ids:
                set_similar_rcsb_ids.update(target_peptide_id2similar_rcsb_ids[target])
            
            if len(set(set_similar_rcsb_ids) & rcsb_train_complex_chain_ids) == 0:
                    novel_target_ids.add(target)
        print(f"Number of new targets: {len(novel_target_ids)}")

        with open(new_target_saving_dir / f"{test_dataset.dataset_name}_new_target_seqs.fasta", "w") as fout:
            for target_id in novel_target_ids:
                fout.write(f">{target_id}\n{target_id2seq[target_id]}\n")

        # exit()
        continue

        # target_id2target_ids = read_alignment(test_dataset.target_to_target_mmseqs_alignment_path)
        # target_peptide_ids = set([line.strip() for line in open(test_dataset.target_peptide_path)])
        # target_peptide_id2similar_rcsb_ids = read_alignment(test_dataset.target_peptide_to_rcsb_train_edit_distance_path, -1)
        

        with open(new_target_saving_dir / f"{test_dataset.dataset_name}_target_seqs.fasta", "w") as fout:
            for target_id in novel_target_ids:
                fout.write(f">{target_id}\n{target_id2seq[target_id]}\n")


if __name__ == "__main__":
    main()
