import json
from collections import Counter, defaultdict

entry_path = "/data/rbg/users/wxsh/esm3/data/antibodies/cov_abdab/ref_seqs/spike_entry_sars1.json"
entry = json.load(open(entry_path))
sequence = entry["sequence"]["value"]
print(len(sequence))


for feat in entry["features"]:
    if feat["type"] == "Region" or feat["type"] == "Domain":
        if feat["description"] == "Receptor-binding domain (RBD)":
            rbd_region = feat
        elif feat["description"] == "BetaCoV S1-NTD":
            ntd_region = feat

rbd_start, rbd_end = rbd_region["location"]["start"]["value"], rbd_region["location"]["end"]["value"]
ntd_start, ntd_end = ntd_region["location"]["start"]["value"], ntd_region["location"]["end"]["value"]
# print(rbd_region)
# print(ntd_region)
# print()
# print(sequence[ntd_start-1:ntd_end])
rbd = sequence[rbd_start-1:rbd_end]
ntd = sequence[ntd_start-1:ntd_end]

with open("/data/rbg/users/wxsh/esm3/data/antibodies/cov_abdab/ref_seqs/SARS-CoV1_WT_S.fasta", "w") as fout:
    fout.write(f">SARS-CoV1\n{sequence}\n")
    fout.write(f">SARS-CoV1_RBD\n{rbd}\n")
    fout.write(f">SARS-CoV1_NTD\n{ntd}\n")

exit()

variants_mutations = defaultdict(list)
# omicron = []
for feat in entry["features"]:
    if feat["type"] == "Natural variant":
        # print(feat["description"])
        if feat["description"].startswith("in strain:"):
            desc = feat["description"].replace("in strain: ", "")
            desc = desc.split(";")[0].strip()
            strains = desc.replace("in strain: ", "")
            strains = [x.strip() for x in strains.split(",")]
            # print(strains)
            # print(feat)
            # location = feat["location"] # ["start"]["value"]
            # alternativeSequence = feat["alternativeSequence"]
            # print(alternativeSequence, location)
            
            for strain in strains:
                variants_mutations[strain].append(feat)
            # break
        else:
            print(feat["description"])
            # print()
        # break
        # variants.append(feat)
        # if "Omicron" in feat["description"] and "BA.1" in feat["description"]:
            # omicron.append(feat)
print(len(variants_mutations))
print(list(variants_mutations.keys()))

# print(variants_mutations["Omicron/BA.1"])
variant_seq = dict()

for variant in variants_mutations:
# for variant in ["Omicron/BA.1", ]:
    feats = variants_mutations[variant]
    print(variant, len(feats))
    seq_list = list(sequence)
    mut_seq_list = list(sequence)
    for feat in feats:
        location = feat["location"] # ["start"]["value"]
        alternativeSequence = feat["alternativeSequence"]
        if len(alternativeSequence) == 0:
            continue
        # print(location)
        if location["start"]["modifier"] == "EXACT":
            start = location["start"]["value"]
        if location["end"]["modifier"] == "EXACT":
            end = location["end"]["value"]
        
        if "originalSequence" not in alternativeSequence:
            print(feat)
        wt = alternativeSequence["originalSequence"]
        mt = "".join(alternativeSequence["alternativeSequences"])
        assert "".join(seq_list[start-1:end]) == wt, f"{start} {end} {seq_list[start-1:end]} {wt}"


        mut_seq_list[start-1]=mt
        for i in range(start, end):
            mut_seq_list[i] = ""
    assert "".join(seq_list) == sequence
    assert "".join(mut_seq_list) != sequence
    # print()
    # print()
    # print(len("".join(mut_seq_list)), len(sequence))
    variant_seq[variant] = ("".join(mut_seq_list), "".join(mut_seq_list[rbd_start-1:rbd_end]), "".join(mut_seq_list[ntd_start-1:ntd_end]))


with open("/data/rbg/users/wxsh/esm3/data/antibodies/cov_abdab/ref_seqs/SARS-CoV1_S_variants.fasta", "w") as fout:
    for v, (full, rbd, ntd) in variant_seq.items():
        fout.write(f">{v}\n{full}\n")
        fout.write(f">{v}_RBD\n{rbd}\n")
        fout.write(f">{v}_NTD\n{ntd}\n")