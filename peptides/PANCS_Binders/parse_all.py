import pandas as pd
from pathlib import Path
from collections import defaultdict, Counter
import re, json

def normalize_target_name(target_name):
    target_name = target_name.replace(" ", "_").replace("(", "_").replace(")", "_").replace(",", "_").replace("\'", "_")
    return target_name

output_dir = Path("processed")
output_dir.mkdir(exist_ok=True, parents=True)

## Read antigens:
antigen_path = "antigens.csv"
antigen_df = pd.read_csv(antigen_path)
antigen_seq = dict()
for target, aa_seq in zip(antigen_df["Target"], antigen_df["AA Sequence"]):
    # target_name = re.sub(r"\s*\([^)]*\)", "", target)
    target = target.strip()
    if target.startswith("SasG") or target.startswith("HSPB1"):
        target_name = target # two segment
    else:
        target_name = re.sub(r"\s*\(\d+-\d+\)", "", target)
    target_name = normalize_target_name(target_name)
    # target_name.replace(" ", "_")
    if target_name in antigen_seq:
        print(target_name, "repeats??")
    else:
        antigen_seq[target_name] = aa_seq

print(f"Read {len(antigen_seq)} antigens.")
with open(output_dir / "pancs_target_seqs.fasta", "w") as fout:
    for antigen_name in antigen_seq:
        fout.write(f">{antigen_name}\n{antigen_seq[antigen_name].strip()}\n")


## Find results for each target
df_target_affibody_affitin_ids = pd.read_csv("target_affibody_affitin_ids.csv")
df_target_affibody_affitin_ids = df_target_affibody_affitin_ids.fillna("")

rename_target = {
    "TRIM": "TRIM21",
    "p65 (structured)": "p65",
    "HlgA toxin": "HlgA",
    "Myc (short)": "Myc (DBD)",
    "Myc (long)": "Myc",
    "OmpF (ec)": "OmpF",
    "HSPB1 (structured)": "HSPB1 (90-171)",
    "HSPB1 (full)": "HSPB1",
    "p53 (long": "p53",
    "p53 (short)": "p53 AD",
    "Calpro het 12aa": "Calprotectin heterodimer with 12aa linker",
    "OmpA (ec)": "OmpA (ec, 84-224)",
    "14-3-3z": "14-3-3z Monomeric",
    "OTUBI": "OTUB1",
    "HIF1a": "HIF1b"    
}

target2affibody_id = dict()
target2affitin_id = dict()
for target_name, affibody_id, affitin_id in zip(df_target_affibody_affitin_ids["Target"], df_target_affibody_affitin_ids["Affibody ID"], df_target_affibody_affitin_ids["Affitin ID"]):
    if target_name in rename_target:
        target_name = rename_target[target_name]

    # target_name = target_name.replace(" ", "_")
    target_name = normalize_target_name(target_name)
    
    if affibody_id:
        if target_name in target2affibody_id:
            print("repeats???", target_name)
        target2affibody_id[target_name] = str(int(affibody_id))
    if affitin_id:
        if target_name in target2affitin_id:
            print("repeats???", target_name)
        target2affitin_id[target_name] = str(int(affitin_id))

    if target_name not in antigen_seq:
        print(f"target without aa seqs: {target_name}")
print(f"Target with affitin: {len(target2affitin_id)}, target with affibody: {len(target2affibody_id)}")


# 
folder_affitin = Path("/data/rbg/users/wxsh/esm3/data/peptides/PANCS_Binders/source_data_2/Affitin Library (Fig 4)")
folder_affibody = Path("/data/rbg/users/wxsh/esm3/data/peptides/PANCS_Binders/source_data_2/Affibody Library (Fig 4)")
read_threshold_low_min = 10
read_threshold_high_min = 100

full_length_affitin = "GSVKVKF???GEEKEVDTSKI??V?R?GK?V?F?YDDNGK?G?G?V?EKDAPKELLDMLARAEREKKL"
full_length_affitin_part =  "???GEEKEVDTSKI??V?R?GK?V?F?YDDNGK?G?G?V?"
full_length_affitin_head = "GSVKVKF"
full_length_affitin_tail = "EKDAPKELLDMLARAEREKKL"

full_length_affibody = "VDNKF?KE???A???I??LPNLN??Q??AF??SL??DPSQSANLLAEAKKLNDAQAPK"
full_length_affibody_part =  "?KE???A???I??LPNLN??Q??AF??SL??"
full_length_affibody_head = "VDNKF"
full_length_affibody_tail = "DPSQSANLLAEAKKLNDAQAPK"

def save_sequences(seqs, output_dir, dataset_name):
    seqs = list(seqs)
    seqs.sort()
    seq2id = dict()
    with open(output_dir / f"{dataset_name}_seqs.fasta", "w") as fout:
        for i, seq in enumerate(seqs):
            fout.write(f">{dataset_name}_{i}\n{seq}\n")
            seq2id[seq] = f"{dataset_name}_{i}"
    return seq2id

def save_records(target2binder_seq2reads, antigen_seq, binder_seq2id, dataset_name):
    records = []
    for target in target2binder_seq2reads:
        for binder_seq in target2binder_seq2reads[target]:
            interactors = []
            interactors.append({"seq": antigen_seq[target], "id": target, "role": "target"})
            interactors.append({"seq": binder_seq, "id": binder_seq2id[binder_seq], "role": "binder"})
            reads, binding = target2binder_seq2reads[target][binder_seq]
            records.append({"interactors": interactors, 
                            "binding": binding, 
                            "attributes": [{"reads": reads}],
                            "assay_id": f"{dataset_name}__{target}"})
    print(f"{dataset_name} records: {len(records)}")
    json.dump(records, open(output_dir / f"{dataset_name}.json", "w"))
    return records

def print_info(target2binder_seq2reads, binder_seqs):
    print("Targets:", len(target2binder_seq2reads))
    print("Binder seqs", len(binder_seqs))
    print("Pairs:", sum([len(x) for x in target2binder_seq2reads.values()]))
    print("Hits (strong):", sum([len([y for y in x.values() if y[1] == 2]) for x in target2binder_seq2reads.values()]))
    print("Hits (weak):", sum([len([y for y in x.values() if y[1] == 1]) for x in target2binder_seq2reads.values()]))


def read_ngs_data(path, full_length_part, full_length_head, full_length_tail, full_length_seq, high_threshold = 0.01, 
                  low_threshold=0.001, quiet=True):
    seq2reads = defaultdict(int)
    number_of_lines = 0
    total_reads = 0
    with open(path) as fin:
        for line in fin:
            number_of_lines += 1
            line = line.strip().split()
            sequence, read, _ = line
            read = int(read)
            seq2reads[sequence] += read
            total_reads += read

    if not quiet:
        print(f"Read {number_of_lines} lines. Number of distinct seqs: {len(seq2reads)}")
        print(f"Total reads {total_reads}.")
        print(Counter([len(seq) for seq in seq2reads]).most_common())
        print(f"Ref length: {len(full_length_seq)}, part length: {len(full_length_part)}, head length: {len(full_length_head)}, tail length: {len(full_length_tail)}")

    ## extract
    regex_pattern = full_length_part.replace('?', '.')
    complete_sequence2reads = defaultdict(int)
    for sequence in seq2reads:
        matches = list(re.finditer(regex_pattern, sequence))
        try:
            assert len(matches) > 0, f"Mismatch: {sequence}, {seq2reads[sequence]}"
            assert len(matches) <= 1, f"More than one substring matches the pattern: {matches}, {sequence}"
            for match in matches:
                start = match.start()
                end = match.end()
                matched_substring = match.group()
                head = sequence[:start]
                tail = sequence[end:]
                if head:
                    assert full_length_head.endswith(head) or head.endswith(full_length_head), f"Head: {head}, full_length_head: {full_length_head}"
                if tail:
                    assert full_length_tail.startswith(tail) or tail.startswith(full_length_tail), f"Tail: {tail}, full_length_tail: {full_length_tail}"
                
                
                complete_sequence = full_length_head + matched_substring + full_length_tail
                assert len(complete_sequence) == len(full_length_seq)
                complete_sequence2reads[complete_sequence] += seq2reads[sequence]
        except Exception as e:
            if not quiet:
                print(e)
    
    read_threshold_high = int(total_reads * high_threshold)
    read_threshold_low = int(total_reads * low_threshold)
    read_threshold_low = max(read_threshold_low, read_threshold_low_min)
    read_threshold_high = max(read_threshold_high, read_threshold_high_min)
    if not quiet:
        print(f"Number of distinct seqs (after matching): {len(complete_sequence2reads)}")
        print(f"read_threshold_low: {read_threshold_low}")
        print(f"read_threshold_high: {read_threshold_high}")
    
    strong_hit_num = 0
    weak_hit_num = 0
    complete_sequence2binding = dict()
    for seq in complete_sequence2reads:
        if complete_sequence2reads[seq] >= read_threshold_high:
            binding = 2
            strong_hit_num += 1
        elif complete_sequence2reads[seq] >= read_threshold_low:
            binding = 1
            weak_hit_num += 1
        else:
            binding = 0
            
        complete_sequence2binding[seq] = (complete_sequence2reads[seq], binding)
        # if complete_sequence2reads[seq] >= threshold:
        #     hit_num += 1
    if not quiet:
        print(f"hits (strong): {strong_hit_num}, hits (weak): {weak_hit_num}")
    return complete_sequence2binding

target2affitin_seq2reads = dict()
affitin_seqs = set()
for target in target2affitin_id:
    # print(">>>", target)
    affitin_path = folder_affitin / (target2affitin_id[target] + ".txt")
    seq2reads = read_ngs_data(
        affitin_path, full_length_affitin_part, full_length_affitin_head, full_length_affitin_tail, full_length_affitin, quiet=True)
    target2affitin_seq2reads[target] = seq2reads
    affitin_seqs.update(seq2reads.keys())

print("Affitin:")
print_info(target2affitin_seq2reads, affitin_seqs)

dataset_name = "pancs_95target_affitin"
affitin_seq2id = save_sequences(affitin_seqs, output_dir, dataset_name)
records_affitin = save_records(target2affitin_seq2reads, antigen_seq, affitin_seq2id, dataset_name)


## Affibody
target2affibody_seq2reads = dict()
affibody_seqs = set()
for target in target2affibody_id:
    if target == "IFNG":
        continue
    # print(">>>", target)
    affibody_path = folder_affibody / (target2affibody_id[target] + ".txt")
    if not affibody_path.exists():
        print(f"FILE DOES NOT EXIST!: {affibody_path}")
        continue
    seq2reads = read_ngs_data(
        affibody_path, full_length_affibody_part, full_length_affibody_head, full_length_affibody_tail, full_length_affibody, quiet=True)
    target2affibody_seq2reads[target] = seq2reads
    affibody_seqs.update(seq2reads.keys())

print("Affibody:")
print_info(target2affibody_seq2reads, affibody_seqs)

dataset_name = "pancs_95target_affibody"
affibody_seq2id = save_sequences(affibody_seqs, output_dir, dataset_name)
records_affibody = save_records(target2affibody_seq2reads, antigen_seq, affibody_seq2id, dataset_name)

### Large library
df_ids_large_lib = pd.read_csv("target_affibody_affitin_ids_large_library.csv")
df_ids_large_lib = df_ids_large_lib.fillna("")

target2affibody_id_large_lib = dict()
target2affitin_id_large_lib = dict()
for target_name, affibody_id, affitin_id in zip(df_ids_large_lib["Target"], df_ids_large_lib["Affibody ID"], df_ids_large_lib["Affitin ID"]):
    if target_name in rename_target:
        target_name = rename_target[target_name]

    # target_name = target_name.replace(" ", "_")
    target_name = normalize_target_name(target_name)
    
    if affibody_id:
        if target_name in target2affibody_id_large_lib:
            print("repeats???", target_name)
        target2affibody_id_large_lib[target_name] = affibody_id
    if affitin_id:
        if target_name in target2affitin_id_large_lib:
            print("repeats???", target_name)
        target2affitin_id_large_lib[target_name] = affitin_id

    if target_name not in antigen_seq:
        print(f"target without aa seqs: {target_name}")

folder_large_library = Path("/data/rbg/users/wxsh/esm3/data/peptides/PANCS_Binders/source_data_2/Large Library PANCS (Fig 5)")

# affibody
target2affibody_seq2reads_large_lib = dict()
affibody_seqs_large_library = set()
for target in target2affibody_id_large_lib:
    affibody_id = target2affibody_id_large_lib[target]
    # print(">>>", target, affibody_id)
    affibody_paths = list((folder_large_library / affibody_id).rglob("Translated_*_GOOD.txt"))
    assert len(affibody_paths) == 1
    for affibody_path in affibody_paths:
        seq2reads = read_ngs_data(
            affibody_path, full_length_affibody_part, full_length_affibody_head, full_length_affibody_tail, full_length_affibody, quiet=True)
        target2affibody_seq2reads_large_lib[target] = seq2reads
        affibody_seqs_large_library.update(seq2reads.keys())

print("Affibody (large library)")
print_info(target2affibody_seq2reads_large_lib, affibody_seqs_large_library)

dataset_name = "pancs_large_lib_affibody"
affibody_seq2id_large_lib = save_sequences(affibody_seqs_large_library, output_dir, dataset_name)
records_affibody_large_lib = save_records(target2affibody_seq2reads_large_lib, antigen_seq, affibody_seq2id_large_lib, dataset_name)

# affitin
target2affitin_seq2reads_large_lib = dict()
affitin_seqs_large_library = set()
for target in target2affitin_id_large_lib: # target2affitin_id:
    affitin_id = target2affitin_id_large_lib[target]
    # print(">>>", target, affitin_id)
    affitin_paths = list((folder_large_library / affitin_id).rglob("Translated_*_GOOD.txt"))
    assert len(affitin_paths) == 1
    for affitin_path in affitin_paths:
        seq2reads = read_ngs_data(
            affitin_path, full_length_affitin_part, full_length_affitin_head, full_length_affitin_tail, full_length_affitin, quiet=True)
        target2affitin_seq2reads_large_lib[target] = seq2reads
        affitin_seqs_large_library.update(seq2reads.keys())

print("Affitin (large library)")
print_info(target2affitin_seq2reads_large_lib, affitin_seqs_large_library)

dataset_name = "pancs_large_lib_affitin"
affitin_seq2id_large_lib = save_sequences(affitin_seqs_large_library, output_dir, dataset_name)
records_affitin_large_lib = save_records(target2affitin_seq2reads_large_lib, antigen_seq, affitin_seq2id_large_lib, dataset_name)

## merge
records = records_affitin + records_affibody + records_affitin_large_lib + records_affibody_large_lib
print(f"Merged records: {len(records)}")
pair2records = defaultdict(list)
for record in records:
    pair2records[(record["interactors"][0]["seq"], record["interactors"][1]["seq"])].append(record)
print(f"Merged records (remove duplicates): {len(pair2records)}")

conflict_count = 0
merged_records = []
for pair in pair2records:
    bindings = set()
    for record in pair2records[pair]:
        binding = record["binding"]
        bindings.add(binding)
    if len(bindings) > 1:
        conflict_count += 1
        continue
    merged_records.extend(pair2records[pair])
print(f"conflict_count: {conflict_count}")
print(f"Merged records (remove conflict, with duplicates): {len(merged_records)}")

with open(output_dir / f"pancs_merged_affibody_affitin_seqs.fasta", "w") as fout:
    for seq in affibody_seq2id:
        fout.write(f">{affibody_seq2id[seq]}\n{seq}\n")
    for seq in affitin_seq2id:
        fout.write(f">{affitin_seq2id[seq]}\n{seq}\n")
    for seq in affibody_seq2id_large_lib:
        fout.write(f">{affibody_seq2id_large_lib[seq]}\n{seq}\n")
    for seq in affitin_seq2id_large_lib:
        fout.write(f">{affitin_seq2id_large_lib[seq]}\n{seq}\n")

print("Assay statistics", Counter([record["assay_id"] for record in merged_records]))
print("Binding statistics", Counter([record["binding"] for record in merged_records]))
print("Targets", len(set([record["interactors"][0]["seq"] for record in records])))
print("Binders", len(set([record["interactors"][1]["seq"] for record in records])))

json.dump(merged_records, open(output_dir / f"pancs_merged.json", "w"))