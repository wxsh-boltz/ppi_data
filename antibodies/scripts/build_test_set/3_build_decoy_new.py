import json
import numpy as np
from collections import defaultdict

def read_cluster(cluster_path):
    cluster_id2seq_ids = defaultdict(list)
    seq_id2cluster_id = dict()
    for line in open(cluster_path):
        cluster_id, seq_id = line.strip().split()
        cluster_id2seq_ids[cluster_id].append(seq_id)
        if seq_id in seq_id2cluster_id:
            assert seq_id2cluster_id[seq_id] == cluster_id
        seq_id2cluster_id[seq_id] = cluster_id
    return cluster_id2seq_ids, seq_id2cluster_id

def read_edit_distance(path):
    id2similar_ids = defaultdict(list)
    with open(path) as fin:
        for line in fin:
            qid, tid, ed = line.strip().split()
            id2similar_ids[qid].append(tid)
    return id2similar_ids


## sabdab:
split = "test"
path = f"/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/sabdab_summary_all_1vh1vl.test_vh_0.7_vl_0.7_an_0.7/data.json"
saving_path = f"/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/sabdab_summary_all_1vh1vl.test_vh_0.7_vl_0.7_an_0.7/data.pos_and_neg.json"
vh_cluster_path = "/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/sabdab_summary_all_1vh1vl.test_vh_0.7_vl_0.7_an_0.7/clusterRes/vh_cluster.tsv"
vl_cluster_path = "/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/sabdab_summary_all_1vh1vl.test_vh_0.7_vl_0.7_an_0.7/clusterRes/vl_cluster.tsv"
antigen_cluster_path = "/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/sabdab_summary_all_1vh1vl.test_vh_0.7_vl_0.7_an_0.7/clusterRes/antigen_cluster.tsv"
antigen_peptide_edit_distance_path = "/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/sabdab_summary_all_1vh1vl.test_vh_0.7_vl_0.7_an_0.7/antigen_peptide_TO_antigen_peptide.m8"

# # ## patent data
# split = "train"
# path = f"/data/rbg/users/wxsh/esm3/data/antibodies/patent/vincent_data_clean_1vh1vl.{split}.pos.json"
# saving_path = f"/data/rbg/users/wxsh/esm3/data/antibodies/patent/vincent_data_clean_1vh1vl.{split}.json"
# vl_cluster_path = "/data/rbg/users/wxsh/esm3/data/antibodies/patent/clusterRes/vl_seqs-min_seq_id_0.7-c_0.25-cov_mode_0_cluster.tsv"
# vh_cluster_path = "/data/rbg/users/wxsh/esm3/data/antibodies/patent/clusterRes/vh_seqs-min_seq_id_0.7-c_0.25-cov_mode_0_cluster.tsv"
# antigen_cluster_path = "/data/rbg/users/wxsh/esm3/data/antibodies/patent/clusterRes/antigen_seqs-min_seq_id_0.7-c_0.25-cov_mode_1_cluster.tsv"

vl_cluster_id2seq_ids, vl_seq_id2cluster_id = read_cluster(vl_cluster_path)
print(f"VL clusters {len(vl_cluster_id2seq_ids)}, seqs {len(vl_seq_id2cluster_id)}")
vh_cluster_id2seq_ids, vh_seq_id2cluster_id = read_cluster(vh_cluster_path)
print(f"VH clusters {len(vh_cluster_id2seq_ids)}, seqs {len(vh_seq_id2cluster_id)}")
antigen_cluster_id2seq_ids, antigen_seq_id2cluster_id = read_cluster(antigen_cluster_path)
print(f"Antigen clusters {len(antigen_cluster_id2seq_ids)}, seqs {len(antigen_seq_id2cluster_id)}")
antigen_peptide2similar_peptides = read_edit_distance(antigen_peptide_edit_distance_path)
print(f"Antigen peptides: {len(antigen_peptide2similar_peptides)}")

def merge_cluster():
    # merge clusters
    merge_cluster_idx = set()
    for p1 in antigen_peptide2similar_peptides:
        for p2 in antigen_peptide2similar_peptides[p1]:
            if antigen_seq_id2cluster_id[p1] != antigen_seq_id2cluster_id[p2]:
                # print(antigen_seq_id2cluster_id[p1], antigen_seq_id2cluster_id[p2])
                cls1 = antigen_seq_id2cluster_id[p1]
                cls2 = antigen_seq_id2cluster_id[p2]
                if cls1 < cls2:
                    merge_cluster_idx.add((cls1, cls2))
                else:
                    merge_cluster_idx.add((cls2, cls1))
            # print(antigen_seq_id2cluster_id[p1] == antigen_seq_id2cluster_id[p2])

    merge_cluster_idx = list(merge_cluster_idx)
    print("Merge clusters: ", merge_cluster_idx)
    for cls_a, cls_b in merge_cluster_idx:
        cls_b_members = antigen_cluster_id2seq_ids.pop(cls_b)
        antigen_cluster_id2seq_ids[cls_a].extend(cls_b_members)
        for seq_b in cls_b_members:
            antigen_seq_id2cluster_id[seq_b] = cls_a
    # print(antigen_cluster_id2seq_ids[merge_cluster_idx[0][0]])
    # print(antigen_seq_id2cluster_id[merge_cluster_idx[0][1]])
    # print(antigen_seq_id2cluster_id[merge_cluster_idx[0][0]])

merge_cluster()

num_neg_samples = 9 # 1 pos, 9 neg
random_seed = 1005
rng = np.random.default_rng(seed=random_seed)

dataset = json.load(open(path))
print(f"Read {len(dataset)} data")


vh_vl_pairs = list() ### TODO: change to LIST please if we want to keep the distributions 
vhh = list()
vll = list() # 
vl_all = list()

vl2seq = dict()
vh2seq = dict()
existing_pos_cluster_paris = set() # used for excluding negative samples
for data in dataset:
    vhs = []
    vls = []
    antigens = []
    for interactor in data["interactors"]:
        if interactor["role"] == "vh":
            vhs.append(interactor["id"])
            vh2seq[interactor["id"]] = interactor["seq"]
        elif interactor["role"] == "vl":
            vls.append(interactor["id"])
            vl2seq[interactor["id"]] = interactor["seq"]
        elif interactor["role"] == "antigen":
            antigens.append(interactor["id"])
    assert len(vhs) <= 1
    assert len(vls) <= 1
    vh = vhs[0] if len(vhs) == 1 else ""
    vl = vls[0] if len(vls) == 1 else ""
    assert len(antigens) > 0
    antigens = tuple(antigens)
    antigen_clusters = tuple([antigen_seq_id2cluster_id[antigen] for antigen in antigens])

    existing_pos_cluster_paris.add((vh_seq_id2cluster_id[vh] if vh else "", vl_seq_id2cluster_id[vl] if vl else "", antigen_clusters))
    for antigen_cluster in antigen_clusters: # could be multiple antigens
        existing_pos_cluster_paris.add((vh_seq_id2cluster_id[vh] if vh else "", vl_seq_id2cluster_id[vl] if vl else "", (antigen_cluster, )))
    if vh and vl:
        existing_pos_cluster_paris.add((vh_seq_id2cluster_id[vh], "", antigen_clusters))
        existing_pos_cluster_paris.add(("", vl_seq_id2cluster_id[vl], antigen_clusters))

    if vh and vl:
        vh_vl_pairs.append((vh, vl))
        vl_all.append(vl)
    elif vh and not vl:
        vhh.append(vh)
    elif vl and not vh:
        vl_all.append(vl)
        vll.append(vl)
print(f'vh_vl_pairs: {len(vh_vl_pairs)}')
print(f'vhh: {len(vhh)}')
print(f'vll: {len(vll)}')
print(f'vl (all): {len(vl_all)}')
print(f"existing_pos_cluster_paris: {len(existing_pos_cluster_paris)}")

# vh_vl_pairs = list(vh_vl_pairs)
vh_vl_pairs.sort()
# vhh = list(vhh)
vhh.sort()
# vl_all = list(vl_all)
vl_all.sort()
vll.sort()

#### Sampling negatives:
def is_in_training_set(neg_vh_cluster, neg_vl_cluster, antigen_clusters, existing_pos_cluster_paris):
    if (neg_vh_cluster, neg_vl_cluster, antigen_clusters) in existing_pos_cluster_paris:
        return True
    for antigen_cluster in antigen_clusters:
        if (neg_vh_cluster, neg_vl_cluster, (antigen_cluster, )) in existing_pos_cluster_paris:
            return True
    return False


train_data_neg = []
existing_neg_paris = set()
for data in dataset:
    data["binding"] = True

    antigens_interactor = []
    antigens = []
    vhs = []
    vls = []
    for interactor in data["interactors"]:
        if interactor["role"] == "antigen":
            antigens_interactor.append(interactor)
            antigens.append(interactor["id"])
        elif interactor["role"] == "vh":
            vhs.append(interactor["id"])
        elif interactor["role"] == "vl":
            vls.append(interactor["id"])
    assert len(vhs) <= 1
    assert len(vls) <= 1
    assert len(antigens_interactor) >= 1
    antigens = tuple(antigens)
    antigen_clusters = tuple([antigen_seq_id2cluster_id[antigen] for antigen in antigens])
    vh = vhs[0] if len(vhs) == 1 else ""
    vl = vls[0] if len(vls) == 1 else ""
    
    neg_vh_vls = set()
    while len(neg_vh_vls) < num_neg_samples:
        if vh and vl:
            neg_vh, neg_vl = vh_vl_pairs[rng.choice(len(vh_vl_pairs), size=1, replace=False)[0]]
            neg_vh_cluster = vh_seq_id2cluster_id[neg_vh]
            neg_vl_cluster = vl_seq_id2cluster_id[neg_vl]

            # if vl_seq_id2cluster_id[vl] != vl_seq_id2cluster_id[neg_vl] and vh_seq_id2cluster_id[vh] != vh_seq_id2cluster_id[neg_vh]:
            if (neg_vh, neg_vl, antigens) not in existing_neg_paris and not is_in_training_set(neg_vh_cluster, neg_vl_cluster, antigen_clusters, existing_pos_cluster_paris):
            # (neg_vh_cluster, neg_vl_cluster, antigen_clusters) not in existing_pos_cluster_paris:
                neg_vh_vls.add((neg_vh, neg_vl))
        elif vh and not vl: # VHH
            neg_vh = vhh[rng.choice(len(vhh), size=1, replace=False)[0]]
            neg_vh_cluster = vh_seq_id2cluster_id[neg_vh]
            neg_vl = ""
            # if vh_seq_id2cluster_id[vh] != vh_seq_id2cluster_id[neg_vh]:
            if (neg_vh, neg_vl, antigens) not in existing_neg_paris and not is_in_training_set(neg_vh_cluster, "", antigen_clusters, existing_pos_cluster_paris):
            # (neg_vh_cluster, "", antigen_clusters) not in existing_pos_cluster_paris:
                neg_vh_vls.add((neg_vh, neg_vl))
        elif vl and not vh: # only light chain, rare cases, not important
            neg_vl = vl_all[rng.choice(len(vl_all), size=1, replace=False)[0]]
            neg_vl_cluster = vl_seq_id2cluster_id[neg_vl]
            neg_vh = ""
            # if vl_seq_id2cluster_id[vl] != vl_seq_id2cluster_id[neg_vl]:
            if (neg_vh, neg_vl, antigens) not in existing_neg_paris and not is_in_training_set("", neg_vl_cluster, antigen_clusters, existing_pos_cluster_paris):
                neg_vh_vls.add((neg_vh, neg_vl))
    
    neg_vh_vls = list(neg_vh_vls)
    neg_vh_vls.sort()
    
    for vh, vl in neg_vh_vls:
        interactors = []
        if vh:
            interactors.append({"id": vh, "seq": vh2seq[vh], "role": "vh"})
        if vl:
            interactors.append({"id": vl, "seq": vl2seq[vl], "role": "vl"})
        interactors.extend(antigens_interactor)
        train_data_neg.append({
            "interactors": interactors, 
            "binding": False, 
            "assay_id": data["assay_id"]
            })
        existing_neg_paris.add((vh, vl, antigens))

print(f"train_data_neg: {len(train_data_neg)}")
print(train_data_neg[0])
# print(train_data_pos[0])
train_data_neg_and_pos = dataset + train_data_neg
print(f"train_data: {len(train_data_neg_and_pos)}")
json.dump(train_data_neg_and_pos, open(saving_path, "w"))


### Check:
neg_cluster_pairs = set()
for data in train_data_neg:
    antigens = []
    vhs = []
    vls = []
    for interactor in data["interactors"]:
        if interactor["role"] == "antigen":
            antigens_interactor.append(interactor)
            antigens.append(interactor["id"])
        elif interactor["role"] == "vh":
            vhs.append(interactor["id"])
        elif interactor["role"] == "vl":
            vls.append(interactor["id"])
    
    antigens = tuple(antigens)
    antigen_clusters = tuple([antigen_seq_id2cluster_id[antigen] for antigen in antigens])
    vh = vhs[0] if len(vhs) == 1 else ""
    vl = vls[0] if len(vls) == 1 else ""

    neg_cluster_pairs.add((vh_seq_id2cluster_id[vh] if vh else "", vl_seq_id2cluster_id[vl] if vl else "", antigen_clusters))
    for antigen_cluster in antigen_clusters: # could be multiple antigens
        neg_cluster_pairs.add((vh_seq_id2cluster_id[vh] if vh else "", vl_seq_id2cluster_id[vl] if vl else "", (antigen_cluster, )))

print(f"neg_cluster_pairs: {len(neg_cluster_pairs)}, overalp with pos: {len(set(neg_cluster_pairs & existing_pos_cluster_paris))}")
# print(len(existing_neg_paris), len(existing_pos_paris), len(existing_neg_paris & existing_pos_paris))
