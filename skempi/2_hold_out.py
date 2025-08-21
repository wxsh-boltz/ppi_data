import pandas as pd
from collections import Counter, defaultdict
from pathlib import Path
from copy import deepcopy
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
import json
from tqdm import tqdm
import numpy as np

def build_holdout_graph(processed_data):
    assay_id2holdout_ids = defaultdict(set)
    for data in processed_data:
        assay_id = data["assay_id"]
        hold_out_proteins = set()
        for attr in data["attributes"]:
            hold_out_proteins.update(attr["interaction_attr"]["Hold_out_proteins"].split(","))
        hold_out_proteins_pdb = {x for x in hold_out_proteins if x not in ("AB/AG", "Pr/PI", "TCR/pMHC")}
        # print(assay_id, hold_out_proteins_pdb)
        assay_id2holdout_ids[assay_id].update(hold_out_proteins_pdb)
        for id in hold_out_proteins_pdb:
            assay_id2holdout_ids[id].add(assay_id)
    return assay_id2holdout_ids

input_path = "processed/skempi_v2_by_pdb.json"
processed_data = json.load(open(input_path))
print(f"Read {len(processed_data)} data")
assay_id2holdout_pdb_assay_ids = build_holdout_graph(processed_data)


hold_out_type_assay_id2records = dict()
for data in processed_data:
    hold_out_types = set()
    for attr in data["attributes"]:
        hold_out_type = attr["interaction_attr"]["Hold_out_type"]
        hold_out_types.add(hold_out_type)
    assert len(hold_out_types) == 1
    hold_out_type = list(hold_out_types)[0]
    if hold_out_type not in hold_out_type_assay_id2records:
        hold_out_type_assay_id2records[hold_out_type] = defaultdict(list)
    hold_out_type_assay_id2records[hold_out_type][data["assay_id"]].append(data)

test = dict()
train = dict()

for hold_out_type in hold_out_type_assay_id2records:
    print(hold_out_type, len(hold_out_type_assay_id2records[hold_out_type]))
    if hold_out_type:
        test.update(hold_out_type_assay_id2records[hold_out_type])
    else:
        train.update(hold_out_type_assay_id2records[hold_out_type])
print("train", len(train), "test", len(test))

initial_test_set = set(test.keys())
holdouts = set()
for assay in initial_test_set:
    holdouts.update(assay_id2holdout_pdb_assay_ids[assay])
holdouts = holdouts - initial_test_set
print(len(holdouts))
while len(holdouts) > 0:
    initial_test_set = holdouts | initial_test_set
    holdouts = set()
    for assay in initial_test_set:
        holdouts.update(assay_id2holdout_pdb_assay_ids[assay])
    holdouts = holdouts - initial_test_set

print(len(holdouts), len(initial_test_set))

for assay in initial_test_set:
    if assay in train:
        test[assay] = train.pop(assay)
print("train_assay", len(train), "test_assay", len(test))
train_record = []
for assay in train:
    train_record.extend(train[assay])

test_record = []
for assay in test:
    test_record.extend(test[assay])

print("train_record", len(train_record), "test_record", len(test_record))
json.dump(train_record, open(input_path.replace(".json", ".train.json"), "w"))
json.dump(test_record, open(input_path.replace(".json", ".test.json"), "w"))

exit()

train_exclude = {}
for assay in test:
    records = test[assay]
    hold_out_proteins = set()
    for record in records:
        for attr in record["attributes"]:
            hold_out_proteins.update(attr["interaction_attr"]["Hold_out_proteins"].split(","))
    
    hold_out_proteins_pdb = {x for x in hold_out_proteins if x not in ("AB/AG", "Pr/PI", "TCR/pMHC")}
    if len(hold_out_proteins_pdb) > 0:
        # print(hold_out_proteins_pdb)
        for id in hold_out_proteins_pdb:
            # assert id not in train, id
            if id in train:
                train_exclude[id] = train.pop(id)

test.update(**train_exclude)
print(len(train), len(test), len(train_exclude))

train_exclude = dict()
for assay in train:
    records = train[assay]
    hold_out_proteins = set()
    for record in records:
        for attr in record["attributes"]:
            hold_out_proteins.update(attr["interaction_attr"]["Hold_out_proteins"].split(","))
    
    # if "AB/AG" in hold_out_proteins_pdb or "Pr/PI" in hold_out_proteins_pdb or "TCR/pMHC" in hold_out_proteins_pdb:
        # print(assay, hold_out_proteins_pdb)
    assert "Pr/PI" not in hold_out_proteins_pdb
    assert "AB/AG" not in hold_out_proteins_pdb
    assert "TCR/pMHC" not in hold_out_proteins_pdb

    hold_out_proteins_pdb = {x for x in hold_out_proteins if x not in ("AB/AG", "Pr/PI", "TCR/pMHC")}
    if len(hold_out_proteins_pdb) > 0:
        for id in hold_out_proteins_pdb:
            # assert id not in train, id
            if id in test:
                train_exclude[assay] = train[assay] # ()

# for data in train:
print(list(train_exclude.keys()))
train = {k: train[k] for k in train if k not in train_exclude}
print(len(train), len(test), len(train_exclude))
test.update(**train_exclude)

### agian?
train_exclude = dict()
for assay in train:
    records = train[assay]
    hold_out_proteins = set()
    for record in records:
        for attr in record["attributes"]:
            hold_out_proteins.update(attr["interaction_attr"]["Hold_out_proteins"].split(","))
    
    # if "AB/AG" in hold_out_proteins_pdb or "Pr/PI" in hold_out_proteins_pdb or "TCR/pMHC" in hold_out_proteins_pdb:
        # print(assay, hold_out_proteins_pdb)
    assert "Pr/PI" not in hold_out_proteins_pdb
    assert "AB/AG" not in hold_out_proteins_pdb
    assert "TCR/pMHC" not in hold_out_proteins_pdb

    hold_out_proteins_pdb = {x for x in hold_out_proteins if x not in ("AB/AG", "Pr/PI", "TCR/pMHC")}
    if len(hold_out_proteins_pdb) > 0:
        for id in hold_out_proteins_pdb:
            # assert id not in train, id
            if id in test:
                train_exclude[assay] = train[assay] # ()

# for data in train:
print(list(train_exclude.keys()))
train = {k: train[k] for k in train if k not in train_exclude}
print(len(train), len(test), len(train_exclude))


test.update(**train_exclude)

### agian?
train_exclude = dict()
for assay in train:
    records = train[assay]
    hold_out_proteins = set()
    for record in records:
        for attr in record["attributes"]:
            hold_out_proteins.update(attr["interaction_attr"]["Hold_out_proteins"].split(","))
    
    # if "AB/AG" in hold_out_proteins_pdb or "Pr/PI" in hold_out_proteins_pdb or "TCR/pMHC" in hold_out_proteins_pdb:
        # print(assay, hold_out_proteins_pdb)
    assert "Pr/PI" not in hold_out_proteins_pdb
    assert "AB/AG" not in hold_out_proteins_pdb
    assert "TCR/pMHC" not in hold_out_proteins_pdb

    hold_out_proteins_pdb = {x for x in hold_out_proteins if x not in ("AB/AG", "Pr/PI", "TCR/pMHC")}
    if len(hold_out_proteins_pdb) > 0:
        for id in hold_out_proteins_pdb:
            # assert id not in train, id
            if id in test:
                train_exclude[assay] = train[assay] # ()

# for data in train:
print(list(train_exclude.keys()))
train = {k: train[k] for k in train if k not in train_exclude}
print(len(train), len(test), len(train_exclude))

exit()
# exit()

path = "skempi_v2.csv"
df = pd.read_csv(path, sep=";")
df = df.fillna("")
# print(df.columns)
print(f"Read {len(df)} data")
print("PDB complex", len(set([x.split("_")[0] for x in df["#Pdb"]])))
print(Counter(df["Hold_out_type"]).most_common())
print(len(set(df["Hold_out_proteins"])))
exit()

counter = []
for hold_out_protein in df["Hold_out_proteins"]:
    counter.extend(hold_out_protein.split(","))
print(len(set(counter)))
print(Counter(counter).most_common())
exit()

all_hold_out_types = set(df["Hold_out_type"])
for hold_out_type in all_hold_out_types:
    print(hold_out_type)
    print(Counter(df[df["Hold_out_type"] == hold_out_type]["Hold_out_proteins"]).most_common())
# for hold_out_type, hold_out_proteins in zip(df["Hold_out_type"], df["Hold_out_proteins"]):
#     print(hold_out_type, hold_out_proteins)
#     exit()


## 1. for each type, keep 30% with more data
## 2. include more based on the Hold_out_proteins
## 3. 