import pandas as pd
import json
from collections import defaultdict
from Bio import SeqIO

hchain_anarci = "processed_old/hchain_seqs.anarci_H.csv"
hchain_anarci_2 = "processed_old/lchain_seqs.anarci_H.csv"
lchain_anarci = "processed_old/lchain_seqs.anarci_KL.csv"
lchain_anarci_2 = "processed_old/hchain_seqs.anarci_KL.csv" # 
input_json_path = "/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/sabdab_summary_all.json"

def extract_results(path):
    df = pd.read_csv(path)

    sequence_cols = list(df.columns)[list(df.columns).index('1'):]
    print(sequence_cols)
    id2seq = defaultdict(list)
    for _, row in df.iterrows():
        row = row.to_dict()
        sequence = "".join([row[pos] for pos in sequence_cols if row[pos] != "-"])
        if not sequence:
            continue
        info = {key : row[key] for key in row if key not in sequence_cols}
        # row["sequence"] = sequence
        id2seq[row["Id"]].append((sequence, info))
        # id2seq[row["Id"]].add(sequence)
        # id2seq[(row["Id"], row["domain_no"])] = (sequence, info)
    print(len(df), len(id2seq))

    # for id in id2seq: # multidomain
    #     if len(id2seq[id]) > 1:
    #         print(id, id2seq[id]) 
    return id2seq

def build_seq2unique_id(id2seq_and_info, saving_path=None):
    seq2ids = defaultdict(list)
    for id in id2seq_and_info:
        seq_and_info_list = id2seq_and_info[id]
        for seq, info in seq_and_info_list:
            seq2ids[seq].append(f"{id}_{info['domain_no']}")
    
    seq2id = dict()
    for s in seq2ids:
        ids = seq2ids[s]
        ids.sort()
        seq2id[s] = ids[0]
    
    if saving_path is not None:
        with open(saving_path, "w") as fout:
            for seq, id in seq2id.items():
                other_ids = seq2ids[seq][1:]
                assert seq2ids[seq][0] == id
                other_ids = "|".join(other_ids)                
                fout.write(f">{id} {other_ids}\n{seq}\n")
    return seq2id

id2vh = extract_results(hchain_anarci)
id2vl = extract_results(lchain_anarci)
id2vl = {**extract_results(lchain_anarci_2), **id2vl}
id2vh = {**extract_results(hchain_anarci_2), **id2vh}

vh_seq2unique_id = build_seq2unique_id(id2vh, "vh_seqs.fasta")
vl_seq2unique_id = build_seq2unique_id(id2vl, "vl_seqs.fasta")
print(f"VH seqs: {len(vh_seq2unique_id)}, VL seqs: {len(vl_seq2unique_id)}")


dataset = json.load(open(input_json_path))
print(len(dataset))
for record in dataset:
    for interactor in record["interactors"]:
        # print(interactor)
        if interactor["role"] == "hchain":
            vhs = id2vh[interactor["id"]]
            domains = [] # 
            for vh in vhs:
                assert (vh[0] in interactor["seq"])
                # print(vh)
                start_idx, end_idx = vh[1]["seqstart_index"], vh[1]["seqend_index"]
                domain_seq = vh[0]
                assert domain_seq == interactor["seq"][start_idx: end_idx + 1]
                domains.append({"seq": domain_seq, "start_idx": start_idx, "end_idx": end_idx, "role": "vh", "id": vh_seq2unique_id[domain_seq]})
            interactor["domains"] = domains     
        elif interactor["role"] == "lchain":
            vls = id2vl[interactor["id"]]
            domains = []
            for vl in vls:
                assert (vl[0] in interactor["seq"])
                start_idx, end_idx = vl[1]["seqstart_index"], vl[1]["seqend_index"]
                domain_seq = vl[0]
                assert domain_seq == interactor["seq"][start_idx: end_idx + 1]
                domains.append({"seq": domain_seq, "start_idx": start_idx, "end_idx": end_idx, "role": "vl", "id": vl_seq2unique_id[domain_seq]})
            interactor["domains"] = domains    
        elif interactor["role"] == "fv":
            domains = []

            vhs = id2vh[interactor["id"]]
            for vh in vhs:
                assert (vh[0] in interactor["seq"])
                start_idx, end_idx = vh[1]["seqstart_index"], vh[1]["seqend_index"]
                domain_seq = vh[0]
                assert domain_seq == interactor["seq"][start_idx: end_idx + 1]
                domains.append({"seq": domain_seq, "start_idx": start_idx, "end_idx": end_idx, "role": "vh", "id": vh_seq2unique_id[domain_seq]})

            vls = id2vl[interactor["id"]]
            for vl in vls:
                assert (vl[0] in interactor["seq"])
                start_idx, end_idx = vl[1]["seqstart_index"], vl[1]["seqend_index"]
                domain_seq = vl[0]
                assert domain_seq == interactor["seq"][start_idx: end_idx + 1]
                domains.append({"seq": domain_seq, "start_idx": start_idx, "end_idx": end_idx, "role": "vl", "id": vl_seq2unique_id[domain_seq]})
            
            # print(domains)
            interactor["domains"] = domains    

            

# merge similar domains
domain_datasets = defaultdict(list)
edge_cases_cnt = 0
for record in dataset:
    # VH(s), VL(s), antigens
    hchains = [chain for chain in record["interactors"] if chain["role"] == "hchain"]
    lchains = [chain for chain in record["interactors"] if chain["role"] == "lchain"]
    antigens = [chain for chain in record["interactors"] if chain["role"] == "antigen"]
    scfv = [chain for chain in record["interactors"] if chain["role"] == "fv"]

    assert len(hchains) <= 1, hchains
    assert len(lchains) <= 1, lchains

    if len(scfv) == 0:
        _is_edge_cases = False
        assert len(hchains) == 1 or len(lchains) == 1 # one heavy chain and one light chain
        hchain_domain_seqs = tuple([x['seq'] for hchain in hchains for x in hchain["domains"]])
        if len(hchain_domain_seqs) > 1:
            # edge_cases_cnt += 1
            _is_edge_cases = True
            # print(hchains)
            # continue
        hchain_domain_seq = hchain_domain_seqs[0] if len(hchain_domain_seqs) == 1 else None
        
        lchain_domain_seqs = tuple([x['seq'] for lchain in lchains for x in lchain["domains"]])
        if len(lchain_domain_seqs) > 1:
            # edge_cases_cnt += 1
            _is_edge_cases = True
            # print(lchains)
            # continue
        lchain_domain_seq = lchain_domain_seqs[0] if len(lchain_domain_seqs) == 1 else None
        
        if _is_edge_cases:
            edge_cases_cnt += 1
            continue


        antigen_seqs = [x["seq"] for x in antigens]
        antigen_seqs.sort()
        domain_datasets[(hchain_domain_seq, lchain_domain_seq, tuple(antigen_seqs))].extend(record["attributes"])
    else:
        vh_domains = []
        vl_domains = []
        assert len(scfv) == 1, "SCFV should be single chain"
        scfv = scfv[0]
        for domain in scfv["domains"]:
            if domain["role"] == "vh":
                vh_domains.append(domain)
            elif domain["role"] == "vl":
                vl_domains.append(domain)
        if len(vh_domains) > 1 or len(vl_domains) > 1:
            edge_cases_cnt += 1
            continue
        
        hchain_domain_seq = vh_domains[0]["seq"] if len(vh_domains) == 1 else None
        lchain_domain_seq = vl_domains[0]["seq"] if len(vl_domains) == 1 else None
        
        antigen_seqs = [x["seq"] for x in antigens]
        antigen_seqs.sort()
        domain_datasets[(hchain_domain_seq, lchain_domain_seq, tuple(antigen_seqs))].extend(record["attributes"])

print(len(domain_datasets))
print("edge_cases_cnt", edge_cases_cnt)

antigen_seq2id = dict()
for record in SeqIO.parse("antigen_seqs.fasta", "fasta"):
    antigen_seq2id[str(record.seq)] = record.id

output_records = []
for (hchain, lchain, antigens), attributes in domain_datasets.items():
    assert hchain != lchain

    interactors = []
    if hchain is not None:
        interactors.append({"seq": hchain, "id": vh_seq2unique_id[hchain], "role": "vh"})
    if lchain is not None:
        interactors.append({"seq": lchain, "id": vl_seq2unique_id[lchain], "role": "vl"})
    for antigen in antigens:
        interactors.append({"seq": antigen, "id": antigen_seq2id[antigen], "role": "antigen"})
    
    output_records.append({"interactors": interactors, "attributes": attributes})

print(len(output_records))

json.dump(output_records, open(input_json_path.replace(".json", "_1vh1vl.json"), "w"))