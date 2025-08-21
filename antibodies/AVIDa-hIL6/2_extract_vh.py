import pandas as pd
import json
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path
import gemmi

saving_dir = Path(".")
hchain_anarci = "AVIDa-hIL6_vhh_seqs.anarci.csv_H.csv"
input_json_path = "AVIDa-hIL6.json"

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
        id2seq[row["Id"].split()[0]].append((sequence, info, cdrs))
    
    # print(len(df), len(id2seq))
    return id2seq

id2vh = parse_anarci_results(hchain_anarci)

def extract_domains(anarci_outputs, interactor, role):
    domains = [] # 
    for seq, info, cdrs in anarci_outputs:
        assert (seq in interactor["seq"]), f"{interactor['id']}, {vh_seq}, {interactor['seq']}"
        start_idx, end_idx = info["seqstart_index"], info["seqend_index"]
        assert seq == interactor["seq"][start_idx: end_idx + 1]
        domains.append({
            "seq": seq, 
            # "seq_pdbx": seq_pdbx,
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
        if interactor["role"] == "vh":
            vhs = id2vh[interactor["id"]]
            interactor["domains"] = extract_domains(vhs, interactor, "vh")

vh_seq2ids = defaultdict(list)
antigen_seq2id = dict()
seq2cdrs = dict()
for record in dataset:
    for interactor in record["interactors"]:
        if interactor['role'] == "antigen":
            antigen_seq2id[interactor["seq"]] = interactor["id"]
            continue
        for domain in interactor["domains"]:
            if domain["role"] == "vh":
                vh_seq2ids[domain["seq"]].append(domain["id"])
                if domain["seq"] not in seq2cdrs:
                    seq2cdrs[domain["seq"]] = domain["cdrs"]
                else:
                    assert seq2cdrs[domain["seq"]] == domain["cdrs"]

print(f"VH seqs: {len(vh_seq2ids)}")
with open(saving_dir / "vh_seqs.fasta", "w") as fout:
    for seq in vh_seq2ids:
        id = vh_seq2ids[seq][0]
        fout.write(f">{id}\n{seq}\n")

# merge similar domains
domain_datasets = defaultdict(list)
fail_to_extract_vhs = []
for record in dataset:
    # VH(s), VL(s), antigens
    hchains = [chain for chain in record["interactors"] if chain["role"] == "vh"]
    antigens = [chain for chain in record["interactors"] if chain["role"] == "antigen"]
    assert len(hchains) == 1, hchains
    hchain_domain_seq = [x['seq'] for hchain in hchains for x in hchain["domains"]] # [0]
    if len(hchain_domain_seq) == 0:
        fail_to_extract_vhs.append(record)
        continue
    assert len(hchain_domain_seq) == 1, hchains
    hchain_domain_seq = hchain_domain_seq[0]
    antigen_seqs = [x["seq"] for x in antigens]
    assert len(antigen_seqs) == 1
    domain_datasets[(hchain_domain_seq, antigen_seqs[0])].append(record)
    

print(f"{len(fail_to_extract_vhs)} records fail to extract VH.")
print(f"Vh dataset: {len(domain_datasets)}")

output_records = []
skip_conflict_binding_label = 0
for (vh_seq, antigen_seq), records in domain_datasets.items():
    interactors = []
    interactors.append({
        "seq": vh_seq,
        "id": vh_seq2ids[vh_seq][0], 
        "role": "vh",
        "cdrs": seq2cdrs[vh_seq]})
    interactors.append({
        "seq": antigen_seq, 
        "id": antigen_seq2id[antigen_seq], 
        "role": "antigen",
        })
    bindings = set(record["binding"] for record in records)
    if len(bindings) > 1:
        print(f"Conflict binding: {bindings}")
        skip_conflict_binding_label += 1
        continue
    
    attributes = []
    for record in records:
        attributes.extend(record["attributes"])
    output_records.append({
        "interactors": interactors, 
        "attributes": attributes,
        "binding": list(bindings)[0],
        "assay_id": records[0]["assay_id"]
        })

print(f"Skip for conflict binding {skip_conflict_binding_label}")
print(f"Records: {len(output_records)}")
json.dump(output_records, open(input_json_path.replace(".json", "_vh.json"), "w"))