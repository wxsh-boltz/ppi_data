import pandas as pd
import json
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path
import gemmi

saving_dir = Path("processed")
hchain_anarci = "processed/hchain_seqs.anarci_H.csv"
hchain_anarci_2 = "processed/lchain_seqs.anarci_H.csv"
lchain_anarci = "processed/lchain_seqs.anarci_KL.csv"
lchain_anarci_2 = "processed/hchain_seqs.anarci_KL.csv" # 
input_json_path = "processed/sabdab_summary_all.json"

def parse_anarci_results(path):
    df = pd.read_csv(path)
    sequence_cols = list(df.columns)[list(df.columns).index('1'):]

    fwr1_cols = sequence_cols[:sequence_cols.index('27')]
    cdr1_cols = sequence_cols[sequence_cols.index('27'):sequence_cols.index('38')+1] # 27-38
    fwr2_cols = sequence_cols[sequence_cols.index('38')+1:sequence_cols.index('56')]
    cdr2_cols = sequence_cols[sequence_cols.index('56'):sequence_cols.index('65')+1] # 56-65
    fwr3_cols = sequence_cols[sequence_cols.index('65')+1:sequence_cols.index('105')] 
    cdr3_cols = sequence_cols[sequence_cols.index('105'):sequence_cols.index('117')+1] # 105-117
    fwr4_cols = sequence_cols[sequence_cols.index('117')+1:]
    
    print("cdr1_cols", cdr1_cols)
    print("cdr2_cols", cdr2_cols)
    print("cdr3_cols", cdr3_cols)

    def get_seqs_from_cols(cols, row):
        return "".join([row[s] for s in cols if row[s] != "-"])

    id2seq = defaultdict(list)
    for _, row in df.iterrows():
        row = row.to_dict()
        sequence = "".join([row[pos] for pos in sequence_cols if row[pos] != "-"])
        if not sequence:
            continue

        cdr1_s, cdr2_s, cdr3_s = get_seqs_from_cols(cdr1_cols, row), get_seqs_from_cols(cdr2_cols, row), get_seqs_from_cols(cdr3_cols, row)
        fwr1_s, fwr2_s, fwr3_s, fwr4_s = get_seqs_from_cols(fwr1_cols, row), get_seqs_from_cols(fwr2_cols, row), get_seqs_from_cols(fwr3_cols, row), get_seqs_from_cols(fwr4_cols, row)

        assert fwr1_s + cdr1_s + fwr2_s + cdr2_s + fwr3_s + cdr3_s + fwr4_s == sequence

        info = {key : row[key] for key in row if key not in sequence_cols}
        cdrs = (cdr1_s, cdr2_s, cdr3_s)
        # domains = {"cdr1": cdr1_s, "cdr2": cdr2_s, "cdr3": cdr3_s, "fwr1": fwr1_s, "fwr2": fwr2_s, "fwr3": fwr3_s, "fwr4": fwr4_s}
        id2seq[row["Id"]].append((sequence, info, cdrs))
    
    # print(len(df), len(id2seq))
    return id2seq

# def build_seq2unique_id(dataset, saving_path=None):
#     vh_pdbx_seq2ids = defaultdict(list)
#     vl_pdbx_seq2ids = defaultdict(list)
#     pdbx_seq2seq = dict()
#     seq2cdrs = dict()
#     for record in dataset:
#         for interactor in record["interactors"]:
#             if "domains" not in interactor:
#                 continue
#             for domain in interactor["domains"]:
#                 if domain["role"] == "vh":
#                     vh_pdbx_seq2ids[domain["seq_pdbx"]].append(domain["id"])
#                     if domain["seq_pdbx"] not in pdbx_seq2seq:
#                         pdbx_seq2seq[domain["seq_pdbx"]] = domain["seq"]
#                     else:
#                         assert pdbx_seq2seq[domain["seq_pdbx"]] == domain["seq"]
#                     if domain["seq"] not in seq2cdrs:
#                         seq2cdrs[domain["seq"]] = domain["cdrs"]
#                     else:
#                         assert seq2cdrs[domain["seq"]] == domain["cdrs"]
#                 elif domain["role"] == "vl":
#                     vl_pdbx_seq2ids[domain["seq_pdbx"]].append(domain["id"])
#                     if domain["seq_pdbx"] not in pdbx_seq2seq:
#                         pdbx_seq2seq[domain["seq_pdbx"]] = domain["seq"]
#                     else:
#                         assert pdbx_seq2seq[domain["seq_pdbx"]] == domain["seq"]
#                     if domain["seq"] not in seq2cdrs:
#                         seq2cdrs[domain["seq"]] = domain["cdrs"]
#                     else:
#                         assert seq2cdrs[domain["seq"]] == domain["cdrs"]

# def build_seq2unique_id(id2seq_and_info, saving_path=None):
#     seq2ids = defaultdict(list)
#     for id in id2seq_and_info:
#         seq_and_info_list = id2seq_and_info[id]
#         for seq, info, domains in seq_and_info_list:
#             seq2ids[seq].append(f"{id}_{info['domain_no']}")
    
#     seq2id = dict()
#     for s in seq2ids:
#         ids = seq2ids[s]
#         ids.sort()
#         seq2id[s] = ids[0]
    
#     if saving_path is not None:
#         with open(saving_path, "w") as fout:
#             for seq, id in seq2id.items():
#                 other_ids = seq2ids[seq][1:]
#                 assert seq2ids[seq][0] == id
#                 other_ids = "|".join(other_ids)                
#                 fout.write(f">{id} {other_ids}\n{seq}\n")
#     return seq2id

id2vh = parse_anarci_results(hchain_anarci)
id2vl = parse_anarci_results(lchain_anarci)
id2vl = {**parse_anarci_results(lchain_anarci_2), **id2vl}
id2vh = {**parse_anarci_results(hchain_anarci_2), **id2vh}

# for id in id2vh:
#     print(id, id2vh[id][-1])
#     break
# for id in id2vl:
#     print(id, id2vl[id][-1])
#     break
# exit()

# vh_seq2unique_id = build_seq2unique_id(id2vh, saving_dir / "vh_seqs.fasta")
# vl_seq2unique_id = build_seq2unique_id(id2vl, saving_dir / "vl_seqs.fasta")
# print(f"VH seqs: {len(vh_seq2unique_id)}, VL seqs: {len(vl_seq2unique_id)}")


def extract_domains(anarci_outputs, interactor, role):
    domains = [] # 
    for seq, info, cdrs in anarci_outputs:
        assert (seq in interactor["seq"]), f"{interactor['id']}, {vh_seq}, {interactor['seq']}"
        start_idx, end_idx = info["seqstart_index"], info["seqend_index"]
        assert seq == interactor["seq"][start_idx: end_idx + 1]
        pdbx_code = gemmi.expand_one_letter_sequence(interactor["seq_pdbx"], gemmi.ResidueKind.AA)
        # print(pdbx_code)
        pdbx_code = pdbx_code[start_idx: end_idx + 1]
        # print(pdbx_code)
        seq_pdbx = gemmi.pdbx_one_letter_code(pdbx_code, gemmi.ResidueKind.AA)
        # print(seq_pdbx)
        domains.append({
            "seq": seq, 
            "seq_pdbx": seq_pdbx,
            # "start_idx": start_idx, 
            # "end_idx": end_idx, 
            "role": role, 
            "id": f"{info['Id']}_{info['domain_no']}", ## here we use the interactor ID pls
            # "id": vh_seq2unique_id[vh_seq], 
            "cdrs": cdrs})
    return domains


dataset = json.load(open(input_json_path))
print(f"Read dataset: {len(dataset)}")
for record in dataset:
    for interactor in record["interactors"]:
        # print(interactor)
        if interactor["role"] == "hchain":
            vhs = id2vh[interactor["id"]]
            interactor["domains"] = extract_domains(vhs, interactor, "vh")
            # print(domains)
            # exit()
            # domains = [] # 
            # for vh_seq, info, cdrs in vhs:
            #     assert (vh_seq.upper() in interactor["seq"].upper()), f"{interactor['id']}, {vh_seq}, {interactor['seq']}, {vh_seq in interactor['seq']}"
            #     start_idx, end_idx = info["seqstart_index"], info["seqend_index"]
            #     assert vh_seq.upper() == interactor["seq"][start_idx: end_idx + 1].upper()

            #     pdbx_code = gemmi.expand_one_letter_sequence(interactor["seq_pdbx"], gemmi.ResidueKind.AA)
            #     print(pdbx_code)
            #     pdbx_code = pdbx_code[start_idx: end_idx + 1]
            #     print(pdbx_code)
            #     seq_pdbx = gemmi.pdbx_one_letter_code(pdbx_code, gemmi.ResidueKind.AA)
            #     print(seq_pdbx)
            #     exit()

            #     domains.append({"seq": vh_seq, "start_idx": start_idx, "end_idx": end_idx, "role": "vh", "id": vh_seq2unique_id[vh_seq], "cdrs": cdrs})
            # interactor["domains"] = domains
        elif interactor["role"] == "lchain":
            vls = id2vl[interactor["id"]]
            interactor["domains"] = extract_domains(vls, interactor, "vl")
            # domains = []
            # for vl_seq, info, cdrs in vls:
            #     assert (vl_seq.upper() in interactor["seq"].upper())
            #     start_idx, end_idx = info["seqstart_index"], info["seqend_index"]
            #     assert vl_seq.upper() == interactor["seq"][start_idx: end_idx + 1].upper()
            #     domains.append({"seq": vl_seq, "start_idx": start_idx, "end_idx": end_idx, "role": "vl", "id": vl_seq2unique_id[vl_seq], "cdrs": cdrs})
            # interactor["domains"] = domains    
        elif interactor["role"] == "fv":
            domains = []

            vhs = id2vh[interactor["id"]]
            domains.extend(extract_domains(vhs, interactor, "vh"))
            
            # for vh_seq, info, cdrs in vhs:
            #     assert (vh_seq.upper() in interactor["seq"].upper())
            #     start_idx, end_idx = info["seqstart_index"], info["seqend_index"]
            #     assert vh_seq.upper() == interactor["seq"][start_idx: end_idx + 1].upper()
            #     domains.append({"seq": vh_seq, "start_idx": start_idx, "end_idx": end_idx, "role": "vh", "id": vh_seq2unique_id[vh_seq], "cdrs": cdrs})

            vls = id2vl[interactor["id"]]
            domains.extend(extract_domains(vls, interactor, "vl"))
            # for vl_seq, info, cdrs in vls:
            #     assert (vl_seq in interactor["seq"])
            #     start_idx, end_idx = info["seqstart_index"], info["seqend_index"]
            #     assert vl_seq == interactor["seq"][start_idx: end_idx + 1]
            #     domains.append({"seq": vl_seq, "start_idx": start_idx, "end_idx": end_idx, "role": "vl", "id": vl_seq2unique_id[vl_seq], "cdrs": cdrs})
            interactor["domains"] = domains    


vh_pdbx_seq2ids = defaultdict(list)
vl_pdbx_seq2ids = defaultdict(list)
antigen_pdbx_seq2id = dict()
pdbx_seq2seq = dict()
seq2cdrs = dict()
for record in dataset:
    for interactor in record["interactors"]:
        if interactor['role'] == "antigen":
            antigen_pdbx_seq2id[interactor["seq_pdbx"]] = interactor["id"]
            pdbx_seq2seq[interactor["seq_pdbx"]] = interactor["seq"]
        # if "domains" not in interactor:
            continue
        for domain in interactor["domains"]:
            if domain["role"] == "vh":
                vh_pdbx_seq2ids[domain["seq_pdbx"]].append(domain["id"])
                if domain["seq_pdbx"] not in pdbx_seq2seq:
                    pdbx_seq2seq[domain["seq_pdbx"]] = domain["seq"]
                else:
                    assert pdbx_seq2seq[domain["seq_pdbx"]] == domain["seq"]
                if domain["seq"] not in seq2cdrs:
                    seq2cdrs[domain["seq"]] = domain["cdrs"]
                else:
                    assert seq2cdrs[domain["seq"]] == domain["cdrs"]
            elif domain["role"] == "vl":
                vl_pdbx_seq2ids[domain["seq_pdbx"]].append(domain["id"])
                if domain["seq_pdbx"] not in pdbx_seq2seq:
                    pdbx_seq2seq[domain["seq_pdbx"]] = domain["seq"]
                else:
                    assert pdbx_seq2seq[domain["seq_pdbx"]] == domain["seq"]
                if domain["seq"] not in seq2cdrs:
                    seq2cdrs[domain["seq"]] = domain["cdrs"]
                else:
                    assert seq2cdrs[domain["seq"]] == domain["cdrs"]

print(f"VH: {len(vh_pdbx_seq2ids)}")
print(f"VL: {len(vl_pdbx_seq2ids)}")

with open(saving_dir / "vh_seqs.fasta", "w") as fout, open(saving_dir / "vh_seqs.pdbx.fasta", "w") as fout_pdbx:
    for pdbx_seq in vh_pdbx_seq2ids:
        seq = pdbx_seq2seq[pdbx_seq]
        id = vh_pdbx_seq2ids[pdbx_seq][0]
        fout_pdbx.write(f">{id}\n{pdbx_seq}\n")
        fout.write(f">{id}\n{seq}\n")
    
with open(saving_dir / "vl_seqs.fasta", "w") as fout, open(saving_dir / "vl_seqs.pdbx.fasta", "w") as fout_pdbx:
    for pdbx_seq in vl_pdbx_seq2ids:
        seq = pdbx_seq2seq[pdbx_seq]
        id = vl_pdbx_seq2ids[pdbx_seq][0]
        fout_pdbx.write(f">{id}\n{pdbx_seq}\n")
        fout.write(f">{id}\n{seq}\n")

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
        hchain_domain_seqs = tuple([x['seq_pdbx'] for hchain in hchains for x in hchain["domains"]])
        if len(hchain_domain_seqs) > 1:
            _is_edge_cases = True
        hchain_domain_seq = hchain_domain_seqs[0] if len(hchain_domain_seqs) == 1 else None
        
        lchain_domain_seqs = tuple([x['seq_pdbx'] for lchain in lchains for x in lchain["domains"]])
        if len(lchain_domain_seqs) > 1:
            _is_edge_cases = True
        lchain_domain_seq = lchain_domain_seqs[0] if len(lchain_domain_seqs) == 1 else None
        
        if _is_edge_cases:
            edge_cases_cnt += 1
            continue

        antigen_seqs = [x["seq_pdbx"] for x in antigens]
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
        
        hchain_domain_seq = vh_domains[0]["seq_pdbx"] if len(vh_domains) == 1 else None
        lchain_domain_seq = vl_domains[0]["seq_pdbx"] if len(vl_domains) == 1 else None
        
        antigen_seqs = [x["seq_pdbx"] for x in antigens]
        antigen_seqs.sort()
        domain_datasets[(hchain_domain_seq, lchain_domain_seq, tuple(antigen_seqs))].extend(record["attributes"])

print(f"Vh/vl dataset: {len(domain_datasets)}")
print("edge_cases_cnt", edge_cases_cnt)

# antigen_pdbx_seq2id = dict()
# for record in SeqIO.parse("antigen_seqs_pdbx.fasta", "fasta"):
#     antigen_pdbx_seq2id[str(record.seq)] = record.id

output_records = []
for (vh_pdbx, vl_pdbx, antigens_pdbx), attributes in domain_datasets.items():
    assert vh_pdbx != vl_pdbx

    interactors = []
    if vh_pdbx is not None:
        vh = pdbx_seq2seq[vh_pdbx]
        interactors.append({
            "seq": vh,
            "id": vh_pdbx_seq2ids[vh_pdbx][0], 
            "role": "vh",
            "seq_pdbx": vh_pdbx, 
            "cdrs": seq2cdrs[vh]})
    if vl_pdbx is not None:
        vl = pdbx_seq2seq[vl_pdbx]
        interactors.append({
            "seq": vl, 
            "id": vl_pdbx_seq2ids[vl_pdbx][0], 
            "role": "vl", 
            "seq_pdbx": vl_pdbx, 
            "cdrs": seq2cdrs[vl]})
    for antigen_pdbx in antigens_pdbx:
        interactors.append({
            "seq": pdbx_seq2seq[antigen_pdbx], 
            "id": antigen_pdbx_seq2id[antigen_pdbx], 
            "role": "antigen",
            "seq_pdbx": antigen_pdbx
            })
    
    output_records.append({"interactors": interactors, "attributes": attributes})

print(len(output_records))

json.dump(output_records, open(input_json_path.replace(".json", "_1vh1vl.json"), "w"))