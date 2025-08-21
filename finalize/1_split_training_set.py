from typing import List
from hydra.utils import instantiate
from omegaconf import OmegaConf
import hydra, json, os
from typing import Any
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
import subprocess
import uuid

MMSEQS_TMP_DIR="tmp"
alnRes_dir = Path("/data/rbg/users/wxsh/esm3/data/finalize/all_alnRes")

def read_antibody_test_set(test_set_path):
    test_set = json.load(open(test_set_path))
    test_vh2complex_id = defaultdict(set)
    test_vl2complex_id = defaultdict(set)
    test_antigen2complex_id = defaultdict(set)
    print(f"Read {len(test_set)} test data.")
    # test_complex = set()
    for i, data in enumerate(test_set):
        # vhs = []
        # vls = []
        # antigens = []
        for interactor in data["interactors"]:
            if interactor["role"] == "vh":
                test_vh2complex_id[interactor["id"]].add(i)
                # vhs.append(interactor["id"])
            elif interactor["role"] == "vl":
                test_vl2complex_id[interactor["id"]].add(i)
                # vls.append(interactor["id"])
            elif interactor["role"] == "antigen":
                test_antigen2complex_id[interactor["id"]].add(i)
                # antigens.append(interactor["id"])
    return test_vh2complex_id, test_vl2complex_id, test_antigen2complex_id
        
        # test_complex.add((tuple(vhs), tuple(vls), tuple(antigens)))
        # for antigen in antigens:
        #     test_complex.add((tuple(vhs), tuple(vls), (antigen, )))

def read_general_ppi_test_set(test_set_path):
    test_set = json.load(open(test_set_path))
    test_seq2complex_id = defaultdict(set)
    print(f"Read {len(test_set)} test data.")
    for i, data in enumerate(test_set):
        for interactor in data["interactors"]:
            test_seq2complex_id[interactor["id"]].add(i)
    return test_seq2complex_id

def read_mini_protein_test_set(test_set_path):
    test_set = json.load(open(test_set_path))
    test_target2complex_id = defaultdict(set)
    test_binder2complex_id = defaultdict(set)
    print(f"Read {len(test_set)} test data.")
    for i, data in enumerate(test_set):
        for interactor in data["interactors"]:
            if interactor["role"] == "target":
                test_target2complex_id[interactor["id"]].add(i)
            if interactor["role"] == "binder":
                test_binder2complex_id[interactor["id"]].add(i)
    return test_target2complex_id, test_binder2complex_id

def read_tcr_pmhc_test_set(test_set_path):
    test_set = json.load(open(test_set_path))
    test_tcrva2complex_id = defaultdict(set)
    test_tcrvb2complex_id = defaultdict(set)
    test_mhc2complex_id = defaultdict(set)
    test_peptide2complex_id = defaultdict(set)
    print(f"Read {len(test_set)} test data.")
    for i, data in enumerate(test_set):
        for interactor in data["interactors"]:
            if interactor["role"] == "mhc":
                test_mhc2complex_id[interactor["id"]].add(i)
            if interactor["role"] == "tcr_va":
                test_tcrva2complex_id[interactor["id"]].add(i)
            if interactor["role"] == "tcr_vb":
                test_tcrvb2complex_id[interactor["id"]].add(i)
            if interactor["role"] == "peptide":
                test_peptide2complex_id[interactor["id"]].add(i)
    return test_mhc2complex_id, test_tcrva2complex_id, test_tcrvb2complex_id, test_peptide2complex_id

EXCLUDE_MODE_1="antigen_and_vh_and_vl"
EXCLUDE_MODE_2="(antigen_and_vh)_or_(antigen_and_vl)"
def exclude_leakage_from_training_set_antibody(
        train_data,
        antigen_id2similar_ids,
        vh_id2similar_ids,
        vl_id2similar_ids,
        test_antigen2complex_id,
        test_vh2complex_id,
        test_vl2complex_id,
        valid_antigen_ids,
        valid_vh_ids,
        valid_vl_ids,
        new_antigen_ids,
        exclude_mode=EXCLUDE_MODE_1
        ):
    assert exclude_mode in (EXCLUDE_MODE_1, EXCLUDE_MODE_2)
    print(f"Read {len(train_data)} train data.")
    new_train_data = []
    skip_for_invalid_antigen = []
    skip_for_invalid_vh = []
    skip_for_invalid_vl = []
    skip_for_new_antigens = []
    for data in train_data:
        vh_similar_complex = set()
        vl_similar_complex = set()
        antigen_similar_complex = set()
        vhs = []
        vls = []
        antigens = []
        antigen_similar_ids = set()
        for interactor in data["interactors"]:
            if interactor["role"] == "antigen":
                ids = antigen_id2similar_ids[interactor["id"]]
                antigen_similar_ids.update(ids)
                for id in ids:
                    antigen_similar_complex.update(test_antigen2complex_id[id])
                antigens.append(interactor["id"])
            elif interactor["role"] == "vh":
                ids = vh_id2similar_ids[interactor["id"]]
                for id in ids:
                    vh_similar_complex.update(test_vh2complex_id[id])
                vhs.append(interactor["id"])
            elif interactor["role"] == "vl":
                vls.append(interactor["id"])
                if vl_id2similar_ids is not None:
                    ids = vl_id2similar_ids[interactor["id"]]
                    for id in ids:
                        vl_similar_complex.update(test_vl2complex_id[id])  
                    
        assert len(antigens) > 0
        assert len(vhs) <= 1
        assert len(vls) <= 1
        assert len(vhs) >= 1 or len(vls) >= 1
        # print(antigens)
        # print(vhs)
        # print(vls)
        # print(test_vh2complex_id[vhs[0]])
        # print(vh_similar_complex)
        # print(vh_id2similar_ids[vhs[0]])
        # # print(vh_id2similar_ids)

        # if len(vh_similar_complex):
        #     print("vh_similar_complex", vh_similar_complex)
        # if len(vl_similar_complex):
        #     print("vl_similar_complex", vl_similar_complex)
        # if len(antigen_similar_complex):
        #     print("antigen_similar_complex", antigen_similar_complex)
        if any(antigen not in valid_antigen_ids for antigen in antigens):
            skip_for_invalid_antigen.append(antigens)
            continue
        if any(vh not in valid_vh_ids for vh in vhs):
            skip_for_invalid_vh.extend(vhs)
            continue 
        if any(vl not in valid_vl_ids for vl in vls):
            skip_for_invalid_vl.extend(vls)
            continue 
        
        if any(similar_antigen in new_antigen_ids for similar_antigen in antigen_similar_ids):
            skip_for_new_antigens.append(antigens)
            continue

        if len(vhs) > 0 and len(vls) > 0:
            if vl_id2similar_ids is not None:
                if exclude_mode == EXCLUDE_MODE_1:
                    if len(vh_similar_complex & vl_similar_complex & antigen_similar_complex) > 0:
                        # print(vh_similar_complex & vl_similar_complex & antigen_similar_complex)
                        continue
                elif exclude_mode == EXCLUDE_MODE_2:
                    if len(antigen_similar_complex & vh_similar_complex) > 0 or len(antigen_similar_complex & vl_similar_complex) > 0:
                        continue
            else: # test set is VHHs
                if len(vh_similar_complex & antigen_similar_complex) > 0:
                    # print(vh_similar_complex & antigen_similar_complex)
                    continue
        elif len(vhs) > 0 and len(vls) == 0: # VHH
            if len(vh_similar_complex & antigen_similar_complex) > 0:
                # print(vh_similar_complex & antigen_similar_complex)
                continue
        elif len(vls) > 0 and len(vhs) == 0: # VL only
            if len(vl_similar_complex & antigen_similar_complex) > 0:
                # print(vl_similar_complex & antigen_similar_complex)
                continue

        # data["split"] = "train"
        new_train_data.append(data)
    print(f"Skip {len(skip_for_invalid_antigen)} interactions for invalid antigens.")
    # print(skip_for_invalid_antigen[:3])
    print(f"Skip {len(skip_for_invalid_vh)} interactions for invalid vhs.")
    print(f"Skip {len(skip_for_invalid_vl)} interactions for invalid vls.")
    print(f"Skip {len(skip_for_new_antigens)} interactions for new antigens.")
    # print(skip_for_new_antigens[:3])
    return new_train_data

def exclude_leakage_from_training_set_tcr_pmhc(
        train_data,
        tcr_va_id2similar_ids,
        tcr_vb_id2similar_ids,
        mhc_id2similar_ids,
        peptide_id2similar_ids,
        test_tcr_va2complex_id,
        test_tcr_vb2complex_id,
        test_mhc2complex_id,
        test_peptide2complex_id,
        valid_mhc_ids,
        valid_tcr_va_ids,
        valid_tcr_vb_ids,
        new_antigen_ids,
        exclude_mode="pair" # or peptide_and_pairs
        ):
    assert exclude_mode in ("pair", "peptide_and_tcr")
    print(f"Read {len(train_data)} train data.")
    new_train_data = []
    skip_for_invalid_seqs = defaultdict(list)
    skip_for_new_antigens = []
    role_train_seq_id2similar_test_ids = {
        "tcr_va": tcr_va_id2similar_ids,
        "tcr_vb": tcr_vb_id2similar_ids,
        "mhc": mhc_id2similar_ids,
        "peptide": peptide_id2similar_ids
    }
    role_test_seq_id2complex_ids = {
        "tcr_va": test_tcr_va2complex_id,
        "tcr_vb": test_tcr_vb2complex_id,
        "mhc": test_mhc2complex_id,
        "peptide": test_peptide2complex_id
    }
    role_valid_ids = {
        "tcr_va": valid_tcr_va_ids,
        "tcr_vb": valid_tcr_vb_ids,
        "mhc": valid_mhc_ids
    }
    for data in train_data:
        role_similar_complex = defaultdict(set)
        role_ids = defaultdict(list)
        role_similar_test_ids = defaultdict(set)

        for interactor in data["interactors"]:
            role = interactor["role"]
            role_ids[role].append(interactor["id"])
            test_ids = role_train_seq_id2similar_test_ids[role][interactor["id"]]
            role_similar_test_ids[role].update(test_ids)
            for id in test_ids:
                role_similar_complex[role].update(role_test_seq_id2complex_ids[role][id])

        for role in ("tcr_va", "tcr_vb", "mhc", "peptide"):
            assert len(role_ids[role]) == 1
        
        _skip = False
        for role in ("tcr_va", "tcr_vb", "mhc"):
            if any(train_id not in role_valid_ids[role] for train_id in role_ids[role]):
                skip_for_invalid_seqs[role].append(role_ids[role])
                _skip = True
                break
        if _skip:
            continue

        _skip = False
        for role in ("tcr_va", "tcr_vb", "mhc", "peptide"):
            if any(test_id in new_antigen_ids for test_id in role_similar_test_ids[role]):
                skip_for_new_antigens.append(role_ids[role])
                _skip = True
                break
        if _skip:
            continue
        
        if exclude_mode == "pair":
            _skip = False
            # print(role_similar_complex)
            # exit()
            flat_role_similar_complex = list(role_similar_complex.values())
            # print(flat_role_similar_complex)
            # print(len(flat_role_similar_complex))
            for i in range(0, len(flat_role_similar_complex)):
                for j in range(i+1, len(flat_role_similar_complex)):
                    if len(flat_role_similar_complex[i] & flat_role_similar_complex[j]) > 0:
                        _skip = True
            if _skip:
                # print(role_similar_complex)
                continue
        elif exclude_mode == "peptide_and_tcr":
            _skip = False
            if len(role_similar_complex["peptide"]) > 0:
                _skip = True
            if len(role_similar_complex["tcr_va"] & role_similar_complex["tcr_vb"]) > 0:
                _skip = True
            if _skip:
                continue
        else:
            raise ValueError(f"Invalid exclude_mode: {exclude_mode}")

        new_train_data.append(data)
    
    for role in skip_for_invalid_seqs:
        print(f"Skip {len(skip_for_invalid_seqs[role])} interactions for invalid {role}.")    
    print(f"Skip {len(skip_for_new_antigens)} interactions for new antigens.")
    return new_train_data

def exclude_leakage_from_training_set_miniproteins(
        train_data,
        target_id2similar_ids,
        binder_id2similar_ids,
        test_target2complex_id,
        test_binder2complex_id,
        valid_target_ids,
        valid_binder_ids,
        new_target_ids,
        new_binder_ids=None # 
        ):
    print(f"Read {len(train_data)} train data.")
    new_train_data = []
    skip_for_invalid_target = []
    skip_for_invalid_binder = []
    skip_for_new_target = []
    skip_for_new_binder = []
    for data in train_data:
        target_similar_complex = set()
        binder_similar_complex = set()
        binders = []
        targets = []
        target_similar_ids = set()
        binder_similar_ids = set()
        for interactor in data["interactors"]:
            if interactor["role"] == "target":
                ids = target_id2similar_ids[interactor["id"]]
                target_similar_ids.update(ids)
                for id in ids:
                    target_similar_complex.update(test_target2complex_id[id])
                targets.append(interactor["id"])
            elif interactor["role"] == "binder":
                ids = binder_id2similar_ids[interactor["id"]]
                binder_similar_ids.update(ids)
                for id in ids:
                    binder_similar_complex.update(test_binder2complex_id[id])
                binders.append(interactor["id"])
                
            
        assert len(binders) == 1
        assert len(targets) == 1

        if any(target not in valid_target_ids for target in targets):
            skip_for_invalid_target.append(targets)
            continue
        if any(binder not in valid_binder_ids for binder in binders):
            skip_for_invalid_binder.extend(binders)
            continue 
        
        if any(similar_target in new_target_ids for similar_target in target_similar_ids):
            skip_for_new_target.append(targets)
            continue

        if new_binder_ids is not None and any(similar_binder in new_binder_ids for similar_binder in binder_similar_ids):
            skip_for_new_binder.append(binders)
            continue

        if len(target_similar_complex & binder_similar_complex) > 0:
            continue

        # data["split"] = "train"
        new_train_data.append(data)
    print(f"Skip {len(skip_for_invalid_target)} interactions for invalid targets.")
    print(f"Skip {len(skip_for_invalid_binder)} interactions for invalid binders.")
    print(f"Skip {len(skip_for_new_target)} interactions for new targets.")
    if new_binder_ids is not None:
        print(f"Skip {len(skip_for_new_binder)} interactions for new binders.")
    return new_train_data

def exclude_leakage_from_training_set_general_ppi(
        train_data,
        train_id2similar_test_ids,
        test_id2test_complex_id,
        valid_seq_ids,
        new_target_ids,
        ):
    print(f"Read {len(train_data)} train data.")
    new_train_data = []
    skip_for_invalid_seq = []
    skip_for_new_target = []
    for data in train_data:
        # assert len(data["interactors"]) == 2
        interactor_similar_complex = [set() for _ in range(len(data["interactors"]))]
        interactor_similar_ids = [set() for _ in range(len(data["interactors"]))]
        for i, interactor in enumerate(data["interactors"]):
            test_ids = train_id2similar_test_ids[interactor["id"]]
            for id in test_ids:
                interactor_similar_complex[i].update(test_id2test_complex_id.get(id, set()))
            interactor_similar_ids[i].update(test_ids)
        
        # print(list(interactor_similar_ids[0])[:2])
        # print(list(interactor_similar_ids[1])[:2])

        if any(interactor["id"] not in valid_seq_ids for interactor in data["interactors"]):
            skip_for_invalid_seq.append(data["interactors"])
            continue
        
        if any(similar_test_id in new_target_ids for similar_test_id in set.union(*interactor_similar_ids)):
            skip_for_new_target.append(data["interactors"])
            continue
        
        #### TODO: this could be problematic. If two interactors are the same, and they are similar to one target, they'll be discarded. 
        _skip = False
        for i in range(len(interactor_similar_complex)):
            for j in range(i+1, len(interactor_similar_complex)):
                if len(interactor_similar_complex[i] & interactor_similar_complex[j]) > 0:
                    _skip = True
        if _skip:
            continue

        # if len(interactor_similar_complex[0] & interactor_similar_complex[1]) > 0:
        #     continue

        new_train_data.append(data)
    print(f"Skip {len(skip_for_invalid_seq)} interactions for invalid seqs.")
    print(f"Skip {len(skip_for_new_target)} interactions for new targets.")
    return new_train_data

def merge_paths(paths):
    if isinstance(paths, str):
        paths = [paths]
    if len(paths) > 1:
        tmp_filename = f"{uuid.uuid4()}.fasta"
        merged_seq_path = Path(MMSEQS_TMP_DIR) / tmp_filename
        with open(merged_seq_path, "w") as fout:
            for path in paths:
                if path is None:
                    continue
                for record in SeqIO.parse(path, "fasta"):
                    fout.write(f">{record.id}\n{record.seq}\n")
        seq_path = merged_seq_path
    else:
        seq_path = paths[0]
        merged_seq_path = None
    return seq_path, merged_seq_path

def run_mmseqs_antibody(
        query_antigen_path, 
        antigen_peptide_path,
        target_antigen_path, 
        query_vh_path, 
        target_vh_path, 
        query_vl_path, 
        target_vl_path, 
        saving_dir, 
        max_peptide_length, 
        force_rerun=False):

    saving_dir = Path(saving_dir)
    max_peptide_length = str(max_peptide_length)

    if not (saving_dir / "antigen_peptide_aln.m8").exists() or force_rerun:
        cmd = [
            "python", "/data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py",
            "--query_fasta_path", antigen_peptide_path,
            "--target_fasta_path", target_antigen_path,
            "--saving_path_alignment", saving_dir / "antigen_peptide_aln.m8",
            "--saving_path_peptides", saving_dir / "train_antigen_peptides.fasta",
            '--max_peptide_length', max_peptide_length,
        ]
        subprocess.run(cmd, check=True)


    if not (saving_dir / "antigen_mode_1.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_antigen_path),
            str(target_antigen_path),
            str(saving_dir / "antigen_mode_1.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "1",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)

    if not (saving_dir / "antigen_mode_2.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_antigen_path),
            str(target_antigen_path),
            str(saving_dir / "antigen_mode_2.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "2",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)

    if not (saving_dir / "vh_mode_0.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_vh_path),
            str(target_vh_path),
            str(saving_dir / "vh_mode_0.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "0",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)

    if not (saving_dir / "vl_mode_0.m8").exists() or force_rerun:
        if query_vl_path is not None and target_vl_path is not None: # nanobody
            cmd = [
                "mmseqs", "easy-search",
                str(query_vl_path),
                str(target_vl_path),
                str(saving_dir / "vl_mode_0.m8"),
                MMSEQS_TMP_DIR,
                "--cov-mode", "0",
                "-c", "0.25",
                "--min-seq-id", "0",
                "--max-seqs", "10000000",
                "-e", "0.1",
                "--spaced-kmer-mode", "0"
            ]
            subprocess.run(cmd, check=True)

def run_mmseqs_miniproteins(
        query_target_path, 
        target_target_path, 
        query_binder_path, 
        target_binder_paths,
        saving_dir, 
        target_peptide_path,
        binder_peptide_path,
        max_peptide_length, 
        force_rerun=False,
        ):
    
    if isinstance(target_binder_paths, str):
        target_binder_paths = [target_binder_paths]

    if len(target_binder_paths) > 1:
        tmp_filename = f"{uuid.uuid4()}.fasta"
        merged_target_binder_path = Path(MMSEQS_TMP_DIR) / tmp_filename
        with open(merged_target_binder_path, "w") as fout:
            for path in target_binder_paths:
                for record in SeqIO.parse(path, "fasta"):
                    fout.write(f">{record.id}\n{record.seq}\n")
        target_binder_path = merged_target_binder_path
    else:
        target_binder_path = target_binder_paths[0]
        merged_target_binder_path = None

    saving_dir = Path(saving_dir)
    max_peptide_length = str(max_peptide_length)

    if not (saving_dir / "target_peptide_aln.m8").exists() or force_rerun:
        cmd = [
            "python", "/data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py",
            "--query_fasta_path", target_peptide_path,
            "--target_fasta_path", target_target_path,
            "--saving_path_alignment", saving_dir / "target_peptide_aln.m8",
            "--saving_path_peptides", saving_dir / "train_target_peptides.fasta",
            '--max_peptide_length', max_peptide_length,
        ]
        subprocess.run(cmd, check=True)
    
    if not (saving_dir / "binder_peptide_aln.m8").exists() or force_rerun:
        cmd = [
            "python", "/data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py",
            "--query_fasta_path", binder_peptide_path,
            "--target_fasta_path", target_binder_path,
            "--saving_path_alignment", saving_dir / "binder_peptide_aln.m8",
            "--saving_path_peptides", saving_dir / "train_binder_peptides.fasta",
            '--max_peptide_length', max_peptide_length,
        ]
        subprocess.run(cmd, check=True)
        
    if not (saving_dir / "target_mode_1.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_target_path),
            str(target_target_path),
            str(saving_dir / "target_mode_1.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "1",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)

    if not (saving_dir / "target_mode_2.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_target_path),
            str(target_target_path),
            str(saving_dir / "target_mode_2.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "2",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)

    if not (saving_dir / "binder_mode_1.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_binder_path),
            str(target_binder_path),
            str(saving_dir / "binder_mode_1.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "1",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)

    if not (saving_dir / "binder_mode_2.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_binder_path),
            str(target_binder_path),
            str(saving_dir / "binder_mode_2.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "2",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)
    
    if merged_target_binder_path is not None:
        os.remove(merged_target_binder_path)

def run_mmseqs_general_ppi(
        query_seq_path,
        target_seq_paths,
        saving_dir, 
        peptide_path,
        max_peptide_length, 
        force_rerun=False,
        ):
    
    if isinstance(target_seq_paths, str):
        target_seq_paths = [target_seq_paths]

    if len(target_seq_paths) > 1:
        tmp_filename = f"{uuid.uuid4()}.fasta"
        merged_seq_path = Path(MMSEQS_TMP_DIR) / tmp_filename
        with open(merged_seq_path, "w") as fout:
            for path in target_seq_paths:
                if path is None:
                    continue
                for record in SeqIO.parse(path, "fasta"):
                    fout.write(f">{record.id}\n{record.seq}\n")
        target_seq_path = merged_seq_path
    else:
        target_seq_path = target_seq_paths[0]
        merged_seq_path = None

    saving_dir = Path(saving_dir)
    max_peptide_length = str(max_peptide_length)

    if not (saving_dir / "peptide_aln.m8").exists() or force_rerun:
        cmd = [
            "python", "/data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py",
            "--query_fasta_path", peptide_path,
            "--target_fasta_path", target_seq_path,
            "--saving_path_alignment", saving_dir / "peptide_aln.m8",
            "--saving_path_peptides", saving_dir / "train_peptides.fasta",
            '--max_peptide_length', max_peptide_length,
        ]
        subprocess.run(cmd, check=True)
    
    if not (saving_dir / "mode_1.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_seq_path),
            str(target_seq_path),
            str(saving_dir / "mode_1.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "1",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)

    if not (saving_dir / "mode_2.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_seq_path),
            str(target_seq_path),
            str(saving_dir / "mode_2.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "2",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)
    
    if merged_seq_path is not None:
        os.remove(merged_seq_path)
    
def run_mmseqs_tcr_pmhc(
        query_tcr_va_path, 
        query_tcr_vb_path,
        query_mhc_path,
        query_peptide_path,
        target_tcr_va_paths,
        target_tcr_vb_paths,
        target_mhc_paths,
        target_peptide_paths,
        saving_dir, 
        max_peptide_length, 
        force_rerun=False):

    saving_dir = Path(saving_dir)
    max_peptide_length = str(max_peptide_length)

    target_tcr_va_path, merged_target_tcr_va_path = merge_paths(target_tcr_va_paths)
    target_tcr_vb_path, merged_target_tcr_vb_path = merge_paths(target_tcr_vb_paths)
    target_mhc_path, merged_target_mhc_path = merge_paths(target_mhc_paths)
    target_peptide_path, merged_target_peptide_path = merge_paths(target_peptide_paths)

    if not (saving_dir / "peptide_aln.m8").exists() or force_rerun:
        cmd = [
            "python", "/data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py",
            "--query_fasta_path", query_peptide_path,
            "--target_fasta_path", target_peptide_path,
            "--saving_path_alignment", saving_dir / "peptide_aln.m8",
            "--saving_path_peptides", saving_dir / "train_peptides.fasta",
            '--max_peptide_length', max_peptide_length,
        ]
        subprocess.run(cmd, check=True)

    if not (saving_dir / "tcr_va_mode_1.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_tcr_va_path),
            str(target_tcr_va_path),
            str(saving_dir / "tcr_va_mode_1.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "1",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)

    if not (saving_dir / "tcr_vb_mode_1.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_tcr_vb_path),
            str(target_tcr_vb_path),
            str(saving_dir / "tcr_vb_mode_1.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "1",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)
    
    if not (saving_dir / "mhc_mode_1.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_mhc_path),
            str(target_mhc_path),
            str(saving_dir / "mhc_mode_1.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "1",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)
    
    if not (saving_dir / "mhc_mode_2.m8").exists() or force_rerun:
        cmd = [
            "mmseqs", "easy-search",
            str(query_mhc_path),
            str(target_mhc_path),
            str(saving_dir / "mhc_mode_2.m8"),
            MMSEQS_TMP_DIR,
            "--cov-mode", "2",
            "-c", "0.25",
            "--min-seq-id", "0",
            "--max-seqs", "10000000",
            "-e", "0.1",
            "--spaced-kmer-mode", "0"
        ]
        subprocess.run(cmd, check=True)

    
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

def read_alignment(
        paths, 
        identity_threshold=0.7, 
        edit_distance_threshold=3, 
        edit_distance_alignment=False,
        skip_self=False):

    query_id2similar_target_ids = defaultdict(list)
    if isinstance(paths, str) or isinstance(paths, Path):
        paths = [paths]
    for path in paths:
        with open(path) as fin:
            for line in fin:
                query_id, target_id, identity = line.strip().split()[:3]
                identity = float(identity)
                if not edit_distance_alignment: # 
                    if identity < identity_threshold:
                        continue
                if edit_distance_alignment:
                    if identity > edit_distance_threshold:
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

def get_peptides(train_data, role=None, max_length=14):
    id2seq = {}
    for data in train_data:
        for interactor in data["interactors"]:
            if role is not None and interactor["role"] != role:
                continue
            if len(interactor["seq"]) <= max_length:
                id2seq[interactor["id"]] = interactor["seq"]
    return id2seq

def read_mmseqs_results_antibody(
        saving_dir, 
        sequence_identity_antigen, 
        sequence_identity_vh, 
        sequence_identity_vl, 
        edit_distance_threshold_antigen):
    antigen2similar_ids = defaultdict(set)
    if (saving_dir / "antigen_mode_1.m8").exists() and (saving_dir / "antigen_mode_2.m8").exists():
        antigen2similar_ids_mode_1 = read_alignment(str(saving_dir / "antigen_mode_1.m8"), sequence_identity_antigen)
        antigen2similar_ids_mode_2 = read_alignment(str(saving_dir / "antigen_mode_2.m8"), sequence_identity_antigen)
        for antigen in antigen2similar_ids_mode_1:
            antigen2similar_ids[antigen].update(antigen2similar_ids_mode_1[antigen])
        for antigen in antigen2similar_ids_mode_2:
            antigen2similar_ids[antigen].update(antigen2similar_ids_mode_2[antigen])
    if (saving_dir / "antigen_peptide_aln.m8").exists():
        antigen2similar_ids_pep = read_alignment(
            saving_dir / "antigen_peptide_aln.m8", 
            edit_distance_alignment=True, 
            identity_threshold=edit_distance_threshold_antigen)
        for antigen in antigen2similar_ids_pep:
            antigen2similar_ids[antigen].update(antigen2similar_ids_pep[antigen])

    if (saving_dir / "vl_mode_0.m8").exists():
        vl2similar_ids = read_alignment(str(saving_dir / "vl_mode_0.m8"), sequence_identity_vh)
    else:
        vl2similar_ids = defaultdict(set)
    if (saving_dir / "vh_mode_0.m8").exists():
        vh2similar_ids = read_alignment(str(saving_dir / "vh_mode_0.m8"), sequence_identity_vl)
    else:
        vh2similar_ids = defaultdict(set)
    return antigen2similar_ids, vh2similar_ids, vl2similar_ids

def read_mmseqs_results_miniproteins(
        saving_dir, 
        sequence_identity_target, 
        sequence_identity_binder, 
        edit_distance_threshold_target,
        edit_distance_threshold_binder):
    target2similar_ids = defaultdict(set)
    if (saving_dir / "target_mode_1.m8").exists() and (saving_dir / "target_mode_2.m8").exists():
        target2similar_ids_mode_1 = read_alignment(str(saving_dir / "target_mode_1.m8"), sequence_identity_target)
        target2similar_ids_mode_2 = read_alignment(str(saving_dir / "target_mode_2.m8"), sequence_identity_target)
        for target in target2similar_ids_mode_1:
            target2similar_ids[target].update(target2similar_ids_mode_1[target])
        for target in target2similar_ids_mode_2:
            target2similar_ids[target].update(target2similar_ids_mode_2[target])
    if (saving_dir / "target_peptide_aln.m8").exists():
        target2similar_ids_pep = read_alignment(
            saving_dir / "target_peptide_aln.m8", 
            edit_distance_alignment=True, 
            identity_threshold=edit_distance_threshold_target)
        for target in target2similar_ids_pep:
            target2similar_ids[target].update(target2similar_ids_pep[target])

    binder2similar_ids = defaultdict(set)
    if (saving_dir / "binder_mode_1.m8").exists() and (saving_dir / "binder_mode_2.m8").exists():
        binder2similar_ids_mode_1 = read_alignment(str(saving_dir / "binder_mode_1.m8"), sequence_identity_binder)
        binder2similar_ids_mode_2 = read_alignment(str(saving_dir / "binder_mode_2.m8"), sequence_identity_binder)
        for binder in binder2similar_ids_mode_1:
            binder2similar_ids[binder].update(binder2similar_ids_mode_1[binder])
        for binder in binder2similar_ids_mode_2:
            binder2similar_ids[binder].update(binder2similar_ids_mode_2[binder])
    if (saving_dir / "binder_peptide_aln.m8").exists():
        binder2similar_ids_pep = read_alignment(
            saving_dir / "binder_peptide_aln.m8", 
            edit_distance_alignment=True, 
            identity_threshold=edit_distance_threshold_binder)
        for binder in binder2similar_ids_pep:
            binder2similar_ids[binder].update(binder2similar_ids_pep[binder])

    return target2similar_ids, binder2similar_ids

def read_mmseqs_results_general_ppi(
        saving_dir, 
        sequence_identity, 
        edit_distance_threshold):
    seq2similar_ids = defaultdict(set)
    if (saving_dir / "mode_1.m8").exists() and (saving_dir / "mode_2.m8").exists():
        similar_ids_mode_1 = read_alignment(str(saving_dir / "mode_1.m8"), sequence_identity)
        similar_ids_mode_2 = read_alignment(str(saving_dir / "mode_2.m8"), sequence_identity)
        for target in similar_ids_mode_1:
            seq2similar_ids[target].update(similar_ids_mode_1[target])
        for target in similar_ids_mode_2:
            seq2similar_ids[target].update(similar_ids_mode_2[target])
    if (saving_dir / "peptide_aln.m8").exists():
        similar_ids_pep = read_alignment(
            saving_dir / "peptide_aln.m8", 
            edit_distance_alignment=True, 
            identity_threshold=edit_distance_threshold)
        for target in similar_ids_pep:
            seq2similar_ids[target].update(similar_ids_pep[target])

    return seq2similar_ids

def read_mmseqs_results_tcr_pmhc(
        saving_dir, 
        sequence_identity_tcr,
        sequence_identity_mhc, 
        edit_distance_threshold_peptide):
    tcr_va_seq2similar_ids = defaultdict(set)
    if (saving_dir / "tcr_va_mode_1.m8").exists():
        similar_ids_mode_1 = read_alignment(str(saving_dir / "tcr_va_mode_1.m8"), sequence_identity_tcr)
        for seq in similar_ids_mode_1:
            tcr_va_seq2similar_ids[seq].update(similar_ids_mode_1[seq])
    
    tcr_vb_seq2similar_ids = defaultdict(set)
    if (saving_dir / "tcr_vb_mode_1.m8").exists():
        similar_ids_mode_1 = read_alignment(str(saving_dir / "tcr_vb_mode_1.m8"), sequence_identity_tcr)
        for seq in similar_ids_mode_1:
            tcr_vb_seq2similar_ids[seq].update(similar_ids_mode_1[seq])

    mhc_seq2similar_ids = defaultdict(set)
    if (saving_dir / "mhc_mode_1.m8").exists() and (saving_dir / "mhc_mode_2.m8").exists():
        similar_ids_mode_1 = read_alignment(str(saving_dir / "mhc_mode_1.m8"), sequence_identity_mhc)
        similar_ids_mode_2 = read_alignment(str(saving_dir / "mhc_mode_2.m8"), sequence_identity_mhc)
        for target in similar_ids_mode_1:
            mhc_seq2similar_ids[target].update(similar_ids_mode_1[target])
        for target in similar_ids_mode_2:
            mhc_seq2similar_ids[target].update(similar_ids_mode_2[target])
    
    peptide_seq2similar_ids = defaultdict(set)
    if (saving_dir / "peptide_aln.m8").exists():
        similar_ids_pep = read_alignment(
            saving_dir / "peptide_aln.m8", 
            edit_distance_alignment=True, 
            identity_threshold=edit_distance_threshold_peptide)
        for target in similar_ids_pep:
            peptide_seq2similar_ids[target].update(similar_ids_pep[target])

    return tcr_va_seq2similar_ids, tcr_vb_seq2similar_ids, mhc_seq2similar_ids, peptide_seq2similar_ids

def save_training_set_antibody(saving_dir, train_data):
    saving_dir = Path(saving_dir)
    saving_dir.mkdir(exist_ok=True, parents=True)
    json.dump(train_data, open(saving_dir / "data.json", "w"))

    role_id2seq = defaultdict(dict)
    for data in train_data:
        for interactor in data["interactors"]:
            if "role" in interactor:
                role_id2seq[interactor["role"]][interactor["id"]] = interactor["seq"]
            else:
                role_id2seq["interactor"][interactor["id"]] = interactor["seq"]
    
    for role in role_id2seq:
        with open(saving_dir / f"{role}_seqs.fasta", "w") as fout:
            for id, seq in role_id2seq[role].items():
                fout.write(f">{id}\n{seq}\n")
    
def statistics(train_data):
    role2ids = defaultdict(set)
    for data in train_data:
        for interactor in data["interactors"]:
            if "role" in interactor:
                role2ids[interactor["role"]].add(interactor["id"])
            else:
                role2ids["interactor"].add(interactor["id"])
    for role in role2ids:
        print(f"{role}: {len(role2ids[role])}")
            
def main():
    rcsb_train_complex_chain_ids = read_rcsb_complexes()
    cfg = OmegaConf.load("datasets.yaml")  # Just load config.yaml
    test_datasets = [instantiate(d) for d in cfg.test_datasets]
    train_datasets = [instantiate(d) for d in cfg.train_datasets]
    
    print(f"Read {len(test_datasets)} test_datasets, and {len(train_datasets)} train_datasets")

    for train_dataset in train_datasets:
        train_data = json.load(open(train_dataset.json_path))
        statistics(train_data)

        saving_dir_train = alnRes_dir / train_dataset.dataset_name
        saving_dir_train.mkdir(exist_ok=True, parents=True)

        print(f"Process training set: {train_dataset.dataset_name} ({train_dataset.dataset_type})")
        
        if train_dataset.dataset_type == "antibody":
            valid_antigen_ids = set()
            valid_antigen_ids.update(read_alignment(train_dataset.alignment_paths.target_to_target).keys())
            # print("valid_antigen_ids:", len(valid_antigen_ids))
            valid_vh_ids = set(read_alignment(train_dataset.alignment_paths.vh_to_vh).keys())
            print("valid_vh_ids:", len(valid_vh_ids))
            valid_vl_ids = set(read_alignment(train_dataset.alignment_paths.vl_to_vl).keys())
            print("valid_vl_ids:", len(valid_vl_ids))

            max_peptide_length = 14
            id2seq =  get_peptides(train_data, "antigen", max_peptide_length)
            valid_antigen_ids.update(id2seq.keys())
            print(f"Number of antigen peptides: {len(id2seq)}")
            antigen_peptide_path = saving_dir_train / "antigen_peptide.fasta"
            with open(antigen_peptide_path, "w") as fout:
                for id, seq in id2seq.items():
                    fout.write(f">{id}\n{seq}\n")
            print("valid_antigen_ids:", len(valid_antigen_ids))

            for test_dataset in test_datasets:
                print(f"\ntrain_dataset: {train_dataset.dataset_name}, test_dataset: {test_dataset.dataset_name}")
                saving_dir = saving_dir_train / test_dataset.dataset_name
                saving_dir.mkdir(exist_ok=True, parents=True)
                
                if test_dataset.dataset_type == "antibody":
                    run_mmseqs_antibody(
                        train_dataset.target_seq_path, 
                        antigen_peptide_path,
                        test_dataset.target_seq_path, 
                        train_dataset.vh_seq_path,
                        test_dataset.vh_seq_path, 
                        train_dataset.vl_seq_path, 
                        test_dataset.vl_seq_path, saving_dir,
                        max_peptide_length
                        )
                    test_vh2complex_id, test_vl2complex_id, test_antigen2complex_id = read_antibody_test_set(test_dataset.json_path)
                    antigen2similar_ids, vh2similar_ids, vl2similar_ids = read_mmseqs_results_antibody(
                        saving_dir, 
                        test_dataset.sequence_identity_threshold.antigen,
                        test_dataset.sequence_identity_threshold.vh,
                        test_dataset.sequence_identity_threshold.vl,
                        test_dataset.edit_distance_threshold.antigen)
                    
                    new_antigen_path = Path(f"new_targets/{test_dataset.dataset_name}_new_target_seqs.fasta")
                    if new_antigen_path.exists():
                        new_antigen_ids = set(record.id for record in SeqIO.parse(new_antigen_path, "fasta"))
                    else:
                        new_antigen_ids = set()
                    
                    new_train_data = exclude_leakage_from_training_set_antibody(
                        train_data, 
                        antigen2similar_ids,
                        vh2similar_ids,
                        vl2similar_ids,
                        test_antigen2complex_id,
                        test_vh2complex_id,
                        test_vl2complex_id,
                        valid_antigen_ids,
                        valid_vh_ids,
                        valid_vl_ids,
                        new_antigen_ids
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "miniprotein":
                    run_mmseqs_antibody(
                        train_dataset.target_seq_path, 
                        antigen_peptide_path,
                        test_dataset.target_seq_path, 
                        train_dataset.vh_seq_path,
                        test_dataset.binder_seq_path, 
                        train_dataset.vl_seq_path, 
                        test_dataset.binder_seq_path, saving_dir,
                        max_peptide_length,
                        )
                    test_target2complex_id, test_binder2complex_id = read_mini_protein_test_set(test_dataset.json_path)
                    antigen2similar_ids, vh2similar_ids, vl2similar_ids = read_mmseqs_results_antibody(
                        saving_dir, 
                        test_dataset.sequence_identity_threshold.target,
                        test_dataset.sequence_identity_threshold.binder,
                        test_dataset.sequence_identity_threshold.binder,
                        test_dataset.edit_distance_threshold.target)
                    new_antigen_path = Path(f"new_targets/{test_dataset.dataset_name}_new_target_seqs.fasta")
                    if new_antigen_path.exists():
                        new_antigen_ids = set(record.id for record in SeqIO.parse(new_antigen_path, "fasta"))
                    else:
                        new_antigen_ids = set()
                    new_train_data = exclude_leakage_from_training_set_antibody(
                        train_data, 
                        antigen2similar_ids,
                        vh2similar_ids,
                        vl2similar_ids,
                        test_target2complex_id,
                        test_binder2complex_id,
                        test_binder2complex_id,
                        valid_antigen_ids,
                        valid_vh_ids,
                        valid_vl_ids,
                        new_antigen_ids,
                        exclude_mode=EXCLUDE_MODE_2
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "general_ppi":
                    run_mmseqs_antibody(
                        train_dataset.target_seq_path, 
                        antigen_peptide_path,
                        test_dataset.seq_path, 
                        train_dataset.vh_seq_path,
                        test_dataset.seq_path, 
                        train_dataset.vl_seq_path, 
                        test_dataset.seq_path, saving_dir,
                        max_peptide_length)

                    test_seq2complex_id = read_general_ppi_test_set(test_dataset.json_path)
                    antigen2similar_ids, vh2similar_ids, vl2similar_ids = read_mmseqs_results_antibody(
                        saving_dir, 
                        test_dataset.sequence_identity_threshold,
                        test_dataset.sequence_identity_threshold,
                        test_dataset.sequence_identity_threshold,
                        test_dataset.edit_distance_threshold)
                    new_train_data = exclude_leakage_from_training_set_antibody(
                        train_data, 
                        antigen2similar_ids,
                        vh2similar_ids,
                        vl2similar_ids,
                        test_seq2complex_id,
                        test_seq2complex_id,
                        test_seq2complex_id,
                        valid_antigen_ids,
                        valid_vh_ids,
                        valid_vl_ids,
                        {},
                        exclude_mode=EXCLUDE_MODE_2
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "tcr_phmc":
                    saving_dir = Path(saving_dir)
                    if not (saving_dir / "antigen_peptide_aln.m8").exists():
                        cmd = [
                            "python", "/data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py",
                            "--query_fasta_path", antigen_peptide_path,
                            "--target_fasta_path", test_dataset.peptide_seq_path,
                            "--saving_path_alignment", saving_dir / "antigen_peptide_aln.m8",
                            "--saving_path_peptides", saving_dir / "train_antigen_peptides.fasta",
                            '--max_peptide_length', "100" # relax, because there will be longer peptides here!
                        ]
                        subprocess.run(cmd, check=True)
                    
                    test_mhc2complex_id, test_tcrva2complex_id, test_tcrvb2complex_id, test_peptide2complex_id = read_tcr_pmhc_test_set(test_dataset.json_path)
                    
                    antigen2similar_ids, vh2similar_ids, vl2similar_ids = read_mmseqs_results_antibody(
                        saving_dir, 
                        None, None, None,
                        # test_dataset.sequence_identity_threshold.target,
                        # test_dataset.sequence_identity_threshold.binder,
                        # test_dataset.sequence_identity_threshold.binder,
                        test_dataset.edit_distance_threshold.peptide)

                    new_antigen_path = Path(test_dataset.peptide_seq_path)
                    if new_antigen_path.exists():
                        new_antigen_ids = set(record.id for record in SeqIO.parse(new_antigen_path, "fasta"))
                    else:
                        new_antigen_ids = set()

                    print(antigen2similar_ids, vh2similar_ids, vl2similar_ids)
                    new_train_data = exclude_leakage_from_training_set_antibody(
                        train_data, 
                        antigen2similar_ids,
                        vh2similar_ids,
                        vl2similar_ids,
                        test_peptide2complex_id,
                        test_tcrva2complex_id,
                        test_tcrvb2complex_id,
                        valid_antigen_ids,
                        valid_vh_ids,
                        valid_vl_ids,
                        new_antigen_ids
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
            
            for data in train_data:
                if "assay_id" not in data:
                    data["assay_id"] = train_dataset.dataset_name
            save_training_set_antibody(train_dataset.json_path.replace(".json", "_train"), train_data)

                
        elif train_dataset.dataset_type == "miniprotein":
            valid_target_ids = set()
            valid_target_ids.update(read_alignment(train_dataset.alignment_paths.target_to_target).keys())
            print("valid_target_ids:", len(valid_target_ids))
            valid_binder_ids = set(read_alignment(train_dataset.alignment_paths.binder_to_binder).keys())
            print("valid_binder_ids:", len(valid_binder_ids))

            max_peptide_length = 14
            target_id2seq =  get_peptides(train_data, "target", max_peptide_length)
            valid_target_ids.update(target_id2seq.keys())
            print(f"Number of target peptides: {len(target_id2seq)}")
            target_peptide_path = saving_dir_train / "target_peptide.fasta"
            with open(target_peptide_path, "w") as fout:
                for id, seq in target_id2seq.items():
                    fout.write(f">{id}\n{seq}\n")
            
            binder_id2seq =  get_peptides(train_data, "binder", max_peptide_length)
            valid_binder_ids.update(binder_id2seq.keys())
            print(f"Number of binder peptides: {len(binder_id2seq)}")
            binder_peptide_path = saving_dir_train / "binder_peptide.fasta"
            with open(binder_peptide_path, "w") as fout:
                for id, seq in binder_id2seq.items():
                    fout.write(f">{id}\n{seq}\n")
            print("valid_binder_ids:", len(valid_binder_ids))

            for test_dataset in test_datasets:
                print(f"\ntrain_dataset: {train_dataset.dataset_name}, test_dataset: {test_dataset.dataset_name}")
                saving_dir = saving_dir_train / test_dataset.dataset_name
                saving_dir.mkdir(exist_ok=True, parents=True)
                
                if test_dataset.dataset_type == "antibody":
                    run_mmseqs_miniproteins(
                        query_target_path=train_dataset.target_seq_path,
                        target_target_path=test_dataset.target_seq_path,
                        query_binder_path=train_dataset.binder_seq_path,
                        target_binder_paths=[test_dataset.vh_seq_path, test_dataset.vh_seq_path],
                        saving_dir=saving_dir,
                        target_peptide_path=target_peptide_path,
                        binder_peptide_path=binder_peptide_path,
                        max_peptide_length=max_peptide_length)
                    test_vh2complex_id, test_vl2complex_id, test_antigen2complex_id = read_antibody_test_set(test_dataset.json_path)
                    target2similar_ids, binder2similar_ids = read_mmseqs_results_miniproteins(
                        saving_dir, 
                        test_dataset.sequence_identity_threshold.antigen,
                        min(test_dataset.sequence_identity_threshold.vh, test_dataset.sequence_identity_threshold.vl),
                        test_dataset.edit_distance_threshold.antigen,
                        test_dataset.edit_distance_threshold.antigen
                        )
                    new_target_path = Path(f"new_targets/{test_dataset.dataset_name}_new_target_seqs.fasta")
                    if new_target_path.exists():
                        new_target_ids = set(record.id for record in SeqIO.parse(new_target_path, "fasta"))
                    else:
                        new_target_ids = set()
                    print("target2similar_ids", len(target2similar_ids))
                    print("binder2similar_ids", len(binder2similar_ids))
                    new_train_data = exclude_leakage_from_training_set_miniproteins(
                        train_data=train_data, 
                        target_id2similar_ids=target2similar_ids,
                        binder_id2similar_ids=binder2similar_ids,
                        test_target2complex_id=test_antigen2complex_id,
                        test_binder2complex_id={**test_vh2complex_id, **test_vl2complex_id},
                        valid_target_ids=valid_target_ids,
                        valid_binder_ids=valid_binder_ids,
                        new_target_ids=new_target_ids
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "miniprotein":
                    run_mmseqs_miniproteins(
                        query_target_path=train_dataset.target_seq_path,
                        target_target_path=test_dataset.target_seq_path,
                        query_binder_path=train_dataset.binder_seq_path,
                        target_binder_paths=test_dataset.binder_seq_path,
                        saving_dir=saving_dir,
                        target_peptide_path=target_peptide_path,
                        binder_peptide_path=binder_peptide_path,
                        max_peptide_length=max_peptide_length)
                    test_target2complex_id, test_binder2complex_id = read_mini_protein_test_set(test_dataset.json_path)
                    target2similar_ids, binder2similar_ids = read_mmseqs_results_miniproteins(
                        saving_dir, 
                        test_dataset.sequence_identity_threshold.target,
                        test_dataset.sequence_identity_threshold.binder,
                        test_dataset.edit_distance_threshold.target,
                        test_dataset.edit_distance_threshold.binder
                        )
                    new_target_path = Path(f"new_targets/{test_dataset.dataset_name}_new_target_seqs.fasta")
                    if new_target_path.exists():
                        new_target_ids = set(record.id for record in SeqIO.parse(new_target_path, "fasta"))
                    else:
                        new_target_ids = set()
                    print("target2similar_ids", len(target2similar_ids))
                    print("binder2similar_ids", len(binder2similar_ids))
                    new_train_data = exclude_leakage_from_training_set_miniproteins(
                        train_data=train_data, 
                        target_id2similar_ids=target2similar_ids,
                        binder_id2similar_ids=binder2similar_ids,
                        test_target2complex_id=test_target2complex_id,
                        test_binder2complex_id=test_binder2complex_id,
                        valid_target_ids=valid_target_ids,
                        valid_binder_ids=valid_binder_ids,
                        new_target_ids=new_target_ids
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "general_ppi":
                    run_mmseqs_miniproteins(
                        query_target_path=train_dataset.target_seq_path,
                        target_target_path=test_dataset.seq_path,
                        query_binder_path=train_dataset.binder_seq_path,
                        target_binder_paths=test_dataset.seq_path,
                        saving_dir=saving_dir,
                        target_peptide_path=target_peptide_path,
                        binder_peptide_path=binder_peptide_path,
                        max_peptide_length=max_peptide_length)
                    test_seq2complex_id = read_general_ppi_test_set(test_dataset.json_path)
                    target2similar_ids, binder2similar_ids = read_mmseqs_results_miniproteins(
                        saving_dir, 
                        test_dataset.sequence_identity_threshold,
                        test_dataset.sequence_identity_threshold,
                        test_dataset.edit_distance_threshold,
                        test_dataset.edit_distance_threshold
                        )
                    print("target2similar_ids", len(target2similar_ids))
                    print("binder2similar_ids", len(binder2similar_ids))
                    new_train_data = exclude_leakage_from_training_set_miniproteins(
                        train_data=train_data,
                        target_id2similar_ids=target2similar_ids,
                        binder_id2similar_ids=binder2similar_ids,
                        test_target2complex_id=test_seq2complex_id,
                        test_binder2complex_id=test_seq2complex_id,
                        valid_target_ids=valid_target_ids,
                        valid_binder_ids=valid_binder_ids,
                        new_target_ids=set()
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "tcr_phmc":
                    saving_dir = Path(saving_dir)
                    if not (saving_dir / "target_peptide_aln.m8").exists():
                        cmd = [
                            "python", "/data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py",
                            "--query_fasta_path", target_peptide_path,
                            "--target_fasta_path", test_dataset.peptide_seq_path,
                            "--saving_path_alignment", saving_dir / "target_peptide_aln.m8",
                            "--saving_path_peptides", saving_dir / "train_target_peptides.fasta",
                            '--max_peptide_length', "100",
                        ]
                        subprocess.run(cmd, check=True)
                    if not (saving_dir / "binder_peptide_aln.m8").exists():
                        cmd = [
                            "python", "/data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py",
                            "--query_fasta_path", binder_peptide_path,
                            "--target_fasta_path", test_dataset.peptide_seq_path,
                            "--saving_path_alignment", saving_dir / "binder_peptide_aln.m8",
                            "--saving_path_peptides", saving_dir / "train_binder_peptides.fasta",
                            '--max_peptide_length', "100",
                        ]
                        subprocess.run(cmd, check=True)
                        
                    test_mhc2complex_id, test_tcrva2complex_id, test_tcrvb2complex_id, test_peptide2complex_id = read_tcr_pmhc_test_set(test_dataset.json_path)
                    target2similar_ids, binder2similar_ids = read_mmseqs_results_miniproteins(
                        saving_dir, 
                        None,
                        None,
                        test_dataset.edit_distance_threshold.peptide,
                        test_dataset.edit_distance_threshold.peptide
                        )

                    new_target_ids = set(record.id for record in SeqIO.parse(test_dataset.peptide_seq_path, "fasta"))
                    new_binder_ids = new_target_ids

                    new_train_data = exclude_leakage_from_training_set_miniproteins(
                        train_data=train_data, 
                        target_id2similar_ids=target2similar_ids,
                        binder_id2similar_ids=binder2similar_ids,
                        test_target2complex_id=test_peptide2complex_id,
                        test_binder2complex_id=test_peptide2complex_id,
                        valid_target_ids=valid_target_ids,
                        valid_binder_ids=valid_binder_ids,
                        new_target_ids=new_target_ids,
                        new_binder_ids=new_binder_ids
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
            
            for data in train_data:
                if "assay_id" not in data:
                    data["assay_id"] = train_dataset.dataset_name
            save_training_set_antibody(train_dataset.json_path.replace(".json", "_train"), train_data)
            
            
        elif train_dataset.dataset_type == "general_ppi":
            valid_seq_ids = set()
            valid_seq_ids.update(read_alignment(train_dataset.alignment_paths.seq_to_seq).keys())
            print("valid_seq_ids:", len(valid_seq_ids))

            max_peptide_length = 14
            peptide_id2seq = get_peptides(train_data, None, max_peptide_length)
            valid_seq_ids.update(peptide_id2seq.keys())
            print(f"Number of peptides: {len(peptide_id2seq)}")
            peptide_path = saving_dir_train / "peptide.fasta"
            with open(peptide_path, "w") as fout:
                for id, seq in peptide_id2seq.items():
                    fout.write(f">{id}\n{seq}\n")
            
            for test_dataset in test_datasets:
                print(f"\ntrain_dataset: {train_dataset.dataset_name}, test_dataset: {test_dataset.dataset_name}")
                saving_dir = saving_dir_train / test_dataset.dataset_name
                saving_dir.mkdir(exist_ok=True, parents=True)
                
                if test_dataset.dataset_type == "antibody":
                    run_mmseqs_general_ppi(
                        query_seq_path=train_dataset.seq_path,
                        target_seq_paths=[test_dataset.target_seq_path, test_dataset.vh_seq_path, test_dataset.vl_seq_path],
                        saving_dir=saving_dir, 
                        peptide_path=peptide_path,
                        max_peptide_length=max_peptide_length)
                    test_vh2complex_id, test_vl2complex_id, test_antigen2complex_id = read_antibody_test_set(test_dataset.json_path)
                    seq2similar_ids = read_mmseqs_results_general_ppi(
                        saving_dir, 
                        min(test_dataset.sequence_identity_threshold.values()),
                        test_dataset.edit_distance_threshold.antigen)                   
                    print("seq2similar_ids", len(seq2similar_ids))
                    new_target_path = Path(f"new_targets/{test_dataset.dataset_name}_new_target_seqs.fasta")
                    if new_target_path.exists():
                        new_target_ids = set(record.id for record in SeqIO.parse(new_target_path, "fasta"))
                    else:
                        new_target_ids = set()
                    new_train_data = exclude_leakage_from_training_set_general_ppi(
                        train_data=train_data,
                        train_id2similar_test_ids=seq2similar_ids,
                        test_id2test_complex_id={**test_vh2complex_id, **test_vl2complex_id, **test_antigen2complex_id},
                        valid_seq_ids=valid_seq_ids,
                        new_target_ids=new_target_ids,
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "miniprotein":
                    run_mmseqs_general_ppi(
                        query_seq_path=train_dataset.seq_path,
                        target_seq_paths=[test_dataset.target_seq_path, test_dataset.binder_seq_path],
                        saving_dir=saving_dir, 
                        peptide_path=peptide_path,
                        max_peptide_length=max_peptide_length)
                    test_target2complex_id, test_binder2complex_id = read_mini_protein_test_set(test_dataset.json_path)
                    # print(list(test_target2complex_id.items())[:3])
                    # print(list(test_binder2complex_id.items())[:3])
                    seq2similar_ids = read_mmseqs_results_general_ppi(
                        saving_dir, 
                        min(test_dataset.sequence_identity_threshold.target, test_dataset.sequence_identity_threshold.binder), 
                        max(test_dataset.edit_distance_threshold.target, test_dataset.edit_distance_threshold.binder))                   
                    # print("seq2similar_ids", len(seq2similar_ids))
                    # print(list(seq2similar_ids.items())[:3])
                    new_target_path = Path(f"new_targets/{test_dataset.dataset_name}_new_target_seqs.fasta")
                    if new_target_path.exists():
                        new_target_ids = set(record.id for record in SeqIO.parse(new_target_path, "fasta"))
                    else:
                        new_target_ids = set()
                    # print(new_target_ids)
                    new_train_data = exclude_leakage_from_training_set_general_ppi(
                        train_data=train_data,
                        train_id2similar_test_ids=seq2similar_ids,
                        test_id2test_complex_id={**test_target2complex_id, **test_binder2complex_id},
                        valid_seq_ids=valid_seq_ids,
                        new_target_ids=new_target_ids,
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "general_ppi":
                    run_mmseqs_general_ppi(
                        query_seq_path=train_dataset.seq_path,
                        target_seq_paths=test_dataset.seq_path,
                        saving_dir=saving_dir, 
                        peptide_path=peptide_path,
                        max_peptide_length=max_peptide_length)
                    test_seq2complex_id = read_general_ppi_test_set(test_dataset.json_path)
                    seq2similar_ids = read_mmseqs_results_general_ppi(
                        saving_dir, 
                        test_dataset.sequence_identity_threshold, 
                        test_dataset.edit_distance_threshold)                   
                    print("seq2similar_ids", len(seq2similar_ids))
                    new_train_data = exclude_leakage_from_training_set_general_ppi(
                        train_data=train_data,
                        train_id2similar_test_ids=seq2similar_ids,
                        test_id2test_complex_id=test_seq2complex_id,
                        valid_seq_ids=valid_seq_ids,
                        new_target_ids=set(),
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "tcr_phmc":
                    saving_dir = Path(saving_dir)
                    if not (saving_dir / "peptide_aln.m8").exists():
                        cmd = [
                            "python", "/data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py",
                            "--query_fasta_path", peptide_path,
                            "--target_fasta_path", test_dataset.peptide_seq_path,
                            "--saving_path_alignment", saving_dir / "peptide_aln.m8",
                            "--saving_path_peptides", saving_dir / "train_peptides.fasta",
                            '--max_peptide_length', "100",
                        ]
                        subprocess.run(cmd, check=True)

                    test_mhc2complex_id, test_tcrva2complex_id, test_tcrvb2complex_id, test_peptide2complex_id = read_tcr_pmhc_test_set(test_dataset.json_path)
                    # print(list(test_mhc2complex_id.items())[:2])
                    # print(list(test_tcrva2complex_id.items())[:2])
                    # print(list(test_tcrvb2complex_id.items())[:2])
                    # print(list(test_peptide2complex_id.items())[:2])
                    seq2similar_ids = read_mmseqs_results_general_ppi(
                        saving_dir, 
                        test_dataset.sequence_identity_threshold, 
                        test_dataset.edit_distance_threshold)

                    new_antigen_path = Path(test_dataset.peptide_seq_path)
                    if new_antigen_path.exists():
                        new_target_ids = set(record.id for record in SeqIO.parse(new_antigen_path, "fasta"))
                    else:
                        new_target_ids = set()

                    print(seq2similar_ids)
                    new_train_data = exclude_leakage_from_training_set_general_ppi(
                        train_data=train_data,
                        train_id2similar_test_ids=seq2similar_ids,
                        test_id2test_complex_id=test_peptide2complex_id,
                        valid_seq_ids=valid_seq_ids,
                        new_target_ids=new_target_ids,
                        )
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
            
            for data in train_data:
                if "assay_id" not in data:
                    data["assay_id"] = train_dataset.dataset_name
            save_training_set_antibody(train_dataset.json_path.replace(".json", "_train"), train_data)
            
        
        elif train_dataset.dataset_type == "tcr_phmc":
            valid_mhc_ids = set(read_alignment(train_dataset.alignment_paths.mhc_to_mhc).keys())
            print("valid_mhc_ids:", len(valid_mhc_ids))
            valid_tcr_va_ids = set(read_alignment(train_dataset.alignment_paths.tcrva_to_tcrva).keys())
            print("valid_tcr_va_ids:", len(valid_tcr_va_ids))
            valid_tcr_vb_ids = set(read_alignment(train_dataset.alignment_paths.tcrvb_to_tcrvb).keys())
            print("valid_tcr_vb_ids:", len(valid_tcr_vb_ids))
            
            
            for test_dataset in test_datasets:
                print(f"\ntrain_dataset: {train_dataset.dataset_name}, test_dataset: {test_dataset.dataset_name}")
                saving_dir = saving_dir_train / test_dataset.dataset_name
                saving_dir.mkdir(exist_ok=True, parents=True)
                
                if test_dataset.dataset_type == "antibody":
                    run_mmseqs_tcr_pmhc(
                        query_tcr_va_path=train_dataset.tcra_seq_path, 
                        query_tcr_vb_path=train_dataset.tcrb_seq_path,
                        query_mhc_path=train_dataset.mhc_seq_path,
                        query_peptide_path=train_dataset.peptide_seq_path,
                        target_tcr_va_paths=[test_dataset.target_seq_path, test_dataset.vh_seq_path, test_dataset.vl_seq_path],
                        target_tcr_vb_paths=[test_dataset.target_seq_path, test_dataset.vh_seq_path, test_dataset.vl_seq_path],
                        target_mhc_paths=[test_dataset.target_seq_path, test_dataset.vh_seq_path, test_dataset.vl_seq_path],
                        target_peptide_paths=[test_dataset.target_seq_path, test_dataset.vh_seq_path, test_dataset.vl_seq_path],
                        saving_dir=saving_dir, 
                        max_peptide_length=100)
                    test_vh2complex_id, test_vl2complex_id, test_antigen2complex_id = read_antibody_test_set(test_dataset.json_path)
                    tcr_va_seq2similar_ids, tcr_vb_seq2similar_ids, mhc_seq2similar_ids, peptide_seq2similar_ids = read_mmseqs_results_tcr_pmhc(
                        saving_dir=saving_dir, 
                        sequence_identity_tcr=min(test_dataset.sequence_identity_threshold.vh, test_dataset.sequence_identity_threshold.vl, test_dataset.sequence_identity_threshold.antigen),
                        sequence_identity_mhc=min(test_dataset.sequence_identity_threshold.vh, test_dataset.sequence_identity_threshold.vl, test_dataset.sequence_identity_threshold.antigen), 
                        edit_distance_threshold_peptide=test_dataset.edit_distance_threshold.antigen
                    )
                    print("tcr_va_seq2similar_ids", len(tcr_va_seq2similar_ids))
                    print("tcr_vb_seq2similar_ids", len(tcr_vb_seq2similar_ids))
                    print("mhc_seq2similar_ids", len(mhc_seq2similar_ids))
                    print("peptide_seq2similar_ids", len(peptide_seq2similar_ids))
                    
                    new_antigen_path = Path(f"new_targets/{test_dataset.dataset_name}_new_target_seqs.fasta")
                    if new_antigen_path.exists():
                        new_antigen_ids = set(record.id for record in SeqIO.parse(new_antigen_path, "fasta"))
                    else:
                        new_antigen_ids = set()
                    
                    new_train_data = exclude_leakage_from_training_set_tcr_pmhc(
                        train_data,
                        tcr_va_seq2similar_ids,
                        tcr_vb_seq2similar_ids,
                        mhc_seq2similar_ids,
                        peptide_seq2similar_ids,
                        {**test_vh2complex_id, **test_vl2complex_id, **test_antigen2complex_id},
                        {**test_vh2complex_id, **test_vl2complex_id, **test_antigen2complex_id},
                        {**test_vh2complex_id, **test_vl2complex_id, **test_antigen2complex_id},
                        {**test_vh2complex_id, **test_vl2complex_id, **test_antigen2complex_id},
                        valid_mhc_ids,
                        valid_tcr_va_ids,
                        valid_tcr_vb_ids,
                        new_antigen_ids,
                        # present_test_roles=("vh", "antigen", ""),
                        exclude_mode="pair")
                        
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "miniprotein":
                    run_mmseqs_tcr_pmhc(
                        query_tcr_va_path=train_dataset.tcra_seq_path, 
                        query_tcr_vb_path=train_dataset.tcrb_seq_path,
                        query_mhc_path=train_dataset.mhc_seq_path,
                        query_peptide_path=train_dataset.peptide_seq_path,
                        target_tcr_va_paths=[test_dataset.target_seq_path, test_dataset.binder_seq_path],
                        target_tcr_vb_paths=[test_dataset.target_seq_path, test_dataset.binder_seq_path],
                        target_mhc_paths=[test_dataset.target_seq_path, test_dataset.binder_seq_path],
                        target_peptide_paths=[test_dataset.target_seq_path, test_dataset.binder_seq_path],
                        saving_dir=saving_dir, 
                        max_peptide_length=100)
                    test_target2complex_id, test_binder2complex_id = read_mini_protein_test_set(test_dataset.json_path)
                    
                    tcr_va_seq2similar_ids, tcr_vb_seq2similar_ids, mhc_seq2similar_ids, peptide_seq2similar_ids = read_mmseqs_results_tcr_pmhc(
                        saving_dir=saving_dir, 
                        sequence_identity_tcr=min(test_dataset.sequence_identity_threshold.target, test_dataset.sequence_identity_threshold.binder),
                        sequence_identity_mhc=min(test_dataset.sequence_identity_threshold.target, test_dataset.sequence_identity_threshold.binder), 
                        edit_distance_threshold_peptide=max(test_dataset.edit_distance_threshold.target, test_dataset.edit_distance_threshold.binder)
                    )
            
                    print("tcr_va_seq2similar_ids", len(tcr_va_seq2similar_ids))
                    print("tcr_vb_seq2similar_ids", len(tcr_vb_seq2similar_ids))
                    print("mhc_seq2similar_ids", len(mhc_seq2similar_ids))
                    print("peptide_seq2similar_ids", len(peptide_seq2similar_ids))
                    
                    new_antigen_path = Path(f"new_targets/{test_dataset.dataset_name}_new_target_seqs.fasta")
                    if new_antigen_path.exists():
                        new_antigen_ids = set(record.id for record in SeqIO.parse(new_antigen_path, "fasta"))
                    else:
                        new_antigen_ids = set()
                    
                    new_train_data = exclude_leakage_from_training_set_tcr_pmhc(
                        train_data,
                        tcr_va_seq2similar_ids,
                        tcr_vb_seq2similar_ids,
                        mhc_seq2similar_ids,
                        peptide_seq2similar_ids,
                        {**test_target2complex_id, **test_binder2complex_id},
                        {**test_target2complex_id, **test_binder2complex_id},
                        {**test_target2complex_id, **test_binder2complex_id},
                        {**test_target2complex_id, **test_binder2complex_id},
                        valid_mhc_ids,
                        valid_tcr_va_ids,
                        valid_tcr_vb_ids,
                        new_antigen_ids,
                        # present_test_roles=("vh", "antigen", ""),
                        exclude_mode="pair")
                        
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "general_ppi":
                    run_mmseqs_tcr_pmhc(
                        query_tcr_va_path=train_dataset.tcra_seq_path, 
                        query_tcr_vb_path=train_dataset.tcrb_seq_path,
                        query_mhc_path=train_dataset.mhc_seq_path,
                        query_peptide_path=train_dataset.peptide_seq_path,
                        target_tcr_va_paths=test_dataset.seq_path,
                        target_tcr_vb_paths=test_dataset.seq_path,
                        target_mhc_paths=test_dataset.seq_path,
                        target_peptide_paths=test_dataset.seq_path,
                        saving_dir=saving_dir, 
                        max_peptide_length=100)
                    test_seq2complex_id = read_general_ppi_test_set(test_dataset.json_path)
                    
                    tcr_va_seq2similar_ids, tcr_vb_seq2similar_ids, mhc_seq2similar_ids, peptide_seq2similar_ids = read_mmseqs_results_tcr_pmhc(
                        saving_dir=saving_dir, 
                        sequence_identity_tcr=test_dataset.sequence_identity_threshold,
                        sequence_identity_mhc=test_dataset.sequence_identity_threshold, 
                        edit_distance_threshold_peptide=test_dataset.edit_distance_threshold
                    )
            
                    print("tcr_va_seq2similar_ids", len(tcr_va_seq2similar_ids))
                    print("tcr_vb_seq2similar_ids", len(tcr_vb_seq2similar_ids))
                    print("mhc_seq2similar_ids", len(mhc_seq2similar_ids))
                    print("peptide_seq2similar_ids", len(peptide_seq2similar_ids))
                    
                    new_train_data = exclude_leakage_from_training_set_tcr_pmhc(
                        train_data,
                        tcr_va_seq2similar_ids,
                        tcr_vb_seq2similar_ids,
                        mhc_seq2similar_ids,
                        peptide_seq2similar_ids,
                        test_seq2complex_id,
                        test_seq2complex_id,
                        test_seq2complex_id,
                        test_seq2complex_id,
                        valid_mhc_ids,
                        valid_tcr_va_ids,
                        valid_tcr_vb_ids,
                        set(), # new_antigen_ids,
                        # present_test_roles=("vh", "antigen", ""),
                        exclude_mode="pair")
                        
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
                elif test_dataset.dataset_type == "tcr_phmc":
                    run_mmseqs_tcr_pmhc(
                        query_tcr_va_path=train_dataset.tcra_seq_path, 
                        query_tcr_vb_path=train_dataset.tcrb_seq_path,
                        query_mhc_path=train_dataset.mhc_seq_path,
                        query_peptide_path=train_dataset.peptide_seq_path,
                        target_tcr_va_paths=test_dataset.tcra_seq_path,
                        target_tcr_vb_paths=test_dataset.tcrb_seq_path,
                        target_mhc_paths=test_dataset.mhc_seq_path,
                        target_peptide_paths=test_dataset.peptide_seq_path,
                        saving_dir=saving_dir, 
                        max_peptide_length=100)

                    test_mhc2complex_id, test_tcrva2complex_id, test_tcrvb2complex_id, test_peptide2complex_id = read_tcr_pmhc_test_set(test_dataset.json_path)

                    tcr_va_seq2similar_ids, tcr_vb_seq2similar_ids, mhc_seq2similar_ids, peptide_seq2similar_ids = read_mmseqs_results_tcr_pmhc(
                        saving_dir=saving_dir, 
                        sequence_identity_tcr=test_dataset.sequence_identity_threshold.tcr,
                        sequence_identity_mhc=test_dataset.sequence_identity_threshold.mhc, 
                        edit_distance_threshold_peptide=test_dataset.edit_distance_threshold.peptide
                    )
            
                    print("tcr_va_seq2similar_ids", len(tcr_va_seq2similar_ids))
                    print("tcr_vb_seq2similar_ids", len(tcr_vb_seq2similar_ids))
                    print("mhc_seq2similar_ids", len(mhc_seq2similar_ids))
                    print("peptide_seq2similar_ids", len(peptide_seq2similar_ids))
                    
                    new_train_data = exclude_leakage_from_training_set_tcr_pmhc(
                        train_data=train_data,
                        tcr_va_id2similar_ids=tcr_va_seq2similar_ids,
                        tcr_vb_id2similar_ids=tcr_vb_seq2similar_ids,
                        mhc_id2similar_ids=mhc_seq2similar_ids,
                        peptide_id2similar_ids=peptide_seq2similar_ids,
                        test_tcr_va2complex_id=test_tcrva2complex_id,
                        test_tcr_vb2complex_id=test_tcrvb2complex_id,
                        test_mhc2complex_id=test_mhc2complex_id,
                        test_peptide2complex_id=test_peptide2complex_id,
                        valid_mhc_ids=valid_mhc_ids,
                        valid_tcr_va_ids=valid_tcr_va_ids,
                        valid_tcr_vb_ids=valid_tcr_vb_ids,
                        new_antigen_ids=set(), # This is included in the exclude-mode!!!
                        # present_test_roles=("vh", "antigen", ""),
                        exclude_mode="peptide_and_tcr")
                        
                    print(f"New training set: {len(new_train_data)}")
                    train_data = new_train_data
            
            for data in train_data:
                if "assay_id" not in data:
                    data["assay_id"] = train_dataset.dataset_name
            save_training_set_antibody(train_dataset.json_path.replace(".json", "_train"), train_data)

            


            
              
            
            
    ### 
    exit()
    for d in test_datasets:
        print(d)
        if d.dataset_type == "antibody":
            dataset = json.load(open(d.json_path))
            print(f"Read {len(dataset)} data")
            target_id2seq = get_target_ids(dataset, d.target_role)
            print(f"Number of targets: {len(target_id2seq)}")
            target_id2similar_rcsb_ids = read_alignment(d.target_to_rcsb_train_mmseqs_alignment_paths)
            target_id2target_ids = read_alignment(d.target_to_target_mmseqs_alignment_path)
            target_pep_ids = set([line.strip() for line in open(d.target_peptide_path)])
            target_pep_id2similar_rcsb_ids = read_alignment(d.target_peptide_to_rcsb_train_edit_distance_path, -1)
            print(len(target_pep_id2similar_rcsb_ids))
            new_target_ids = []
            for target_id in target_id2seq:
                overlap_with_rcsb = rcsb_train_complex_chain_ids & set(target_id2similar_rcsb_ids[target_id])
                if target_id in target_pep_ids:
                    overlap_with_rcsb = overlap_with_rcsb | set(target_pep_id2similar_rcsb_ids[target_id])
                if (len(target_id2target_ids[target_id]) > 0 or target_id in target_pep_ids) and len(overlap_with_rcsb) == 0:
                    new_target_ids.append(target_id)
                
            print(f"New targets: {len(new_target_ids)}")
            
            #### read alignment to the rcsb_train




if __name__ == "__main__":
    main()
