import pandas as pd
from collections import Counter, defaultdict
from pathlib import Path
from copy import deepcopy
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
import json
from tqdm import tqdm
import gemmi
from foldeverything.data import const

def residue_to_char(residue):
    if is_aa(residue):  
        resname = residue.get_resname()
        aa_name_boltz_extended = const.protein_letters_3to1_extended.get(resname, "X")
        if aa_name_boltz_extended == "X":
            aa_name = gemmi.find_tabulated_residue(resname).one_letter_code
        else:
            aa_name = aa_name_boltz_extended
        if not aa_name.strip():
            return "X"
        else:
            return aa_name      
        # aa_name_boltz_extended = const.protein_letters_3to1_extended.get(resname, "X")
        # if aa_name_boltz_extended == "X":
        #     aa_name = gemmi.find_tabulated_residue(aa).one_letter_code
        # else:
        #     aa_name = aa_name_boltz_extended
        # if not aa_name.strip():
        #     return "X"
        # else:
        #     return aa_name
        # # return protein_letters_3to1[resname]
        # # except KeyError:
        #     # print(resname)
        #     # return 'X'  # unknown standard AA
    # elif is_aa(residue, standard=False):
    #     print("not aa?", residue)
    #     return "X" # transfer non-standard AA to X
    return ""  # skip non-AA residues

def extract_residues_from_pdb(pdb_file_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file_path)
    chain_residues = {}
    for model in structure:
        for chain in model:
            residues = []
            for residue in chain:
                aa = residue_to_char(residue)
                residues.append(aa)
            chain_residues[chain.id] = residues

    return chain_residues

path = "skempi_v2.csv"
pdb_dir = "PDBs"
output_dir = Path("processed")
output_dir.mkdir(exist_ok=True, parents=True)

df = pd.read_csv(path, sep=";")
df = df.fillna("")
print(f"Read {len(df)} data")
print("PDB complex", len(set([x.split("_")[0] for x in df["#Pdb"]])))
print(Counter(df["Hold_out_type"]).most_common())

def parse_pdb(folder):
    # pdb_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
    pdb_files = list(folder.glob('*.pdb'))
    all_sequences = dict()
    for pdb_file in tqdm(pdb_files):
        # pdb_path = folder / pdb_file
        pdb_id = pdb_file.stem
        if pdb_id.startswith("."):
            continue
        try:
            sequences = extract_residues_from_pdb(pdb_file)
            all_sequences[pdb_id] = sequences
        except Exception as e:
            print(pdb_id, e)
    return all_sequences

all_sequences_wt = parse_pdb(Path(pdb_dir))
print(f"Parse PDB:", len(all_sequences_wt))
# print(all_sequences_wt)
with open(output_dir / "wt_unique_sequences.fasta", "w") as fout:
    for pdb_id in all_sequences_wt:
        for chain_id in all_sequences_wt[pdb_id]:
            seq = "".join(all_sequences_wt[pdb_id][chain_id])
            if seq:
                fout.write(f">{pdb_id}_{chain_id}\n{seq}\n")

for pdb, mutations in zip(df["#Pdb"], df["Mutation(s)_cleaned"]):
    pdb_id, chains_a, chains_b = pdb.split("_")
    mutations = mutations.split(",")
    for mut in mutations:
        wt = mut[0]
        chain_id = mut[1]
        pos = int(mut[2:-1])
        mt = mut[-1]
        chain_seq = all_sequences_wt[pdb_id][chain_id]
        assert chain_seq[pos-1] == wt
    # if pdb == "4NKQ_C_AB":
    #     print(all_sequences_wt[pdb_id]['B'])
        # if chain_seq[pos-1] != wt:
        #     print(all_sequences_wt[pdb_id])
        #     print(pdb, mut)
        #     print(chain_seq, chain_seq[pos-1])
        #     exit()

pdb2records = dict()
seq2ids = defaultdict(list)

df_wt = df.drop(columns=[
    "Mutation(s)_PDB", "Mutation(s)_cleaned", 
    "iMutation_Location(s)", "Affinity_mut (M)",
    "Affinity_mut_parsed", 'kon_mut (M^(-1)s^(-1))', 'kon_mut_parsed', 'koff_mut (s^(-1))', 'koff_mut_parsed',
    'dH_mut (kcal mol^(-1))', 'dS_mut (cal mol^(-1) K^(-1))', 'dS_wt (cal mol^(-1) K^(-1))'
])
for index, row in df_wt.iterrows():
    row_dict = row.to_dict()
    pdb = row_dict["#Pdb"]
    if pdb not in pdb2records:
        pdb2records[pdb] = defaultdict(list)
    pdb_id, chains_1, chains_2 = pdb.split("_") 
    
    interactor_1 = [] # could have multiple chains
    interactor_2 = []
    for chain_id in chains_1:
        seq = "".join(all_sequences_wt[pdb_id][chain_id])
        id = f"{pdb_id}_{chain_id}"
        interactor_1.append({
            "pdb_chain": f"{pdb_id}_{chain_id}",
            "mutations": [],
            "seq": seq,
        })
        seq2ids[seq].append(id)
    
    for chain_id in chains_2:
        seq = "".join(all_sequences_wt[pdb_id][chain_id])
        id = f"{pdb_id}_{chain_id}"
        interactor_2.append({
            "pdb_chain": f"{pdb_id}_{chain_id}",
            "mutations": [],
            "seq": seq,
        }) 
        seq2ids[seq].append(id)
    
    interactor_1 = sorted(interactor_1, key=lambda x: x["seq"])
    interactor_2 = sorted(interactor_2, key=lambda x: x["seq"])

    if tuple([x["seq"] for x in interactor_1]) > tuple([x["seq"] for x in interactor_2]):
        interactor_1, interactor_2 = interactor_2, interactor_1
        
    pdb2records[pdb][(tuple(x["seq"] for x in interactor_1), tuple(x["seq"] for x in interactor_2))].append(
        {
            "interactors": [interactor_1, interactor_2], 
            "interaction_attr": row_dict,
            "wt": True
        })

for index, row in df.iterrows():
    row_dict = row.to_dict()
    pdb = row_dict["#Pdb"]
    if pdb not in pdb2records:
        pdb2records[pdb] = defaultdict(list)
    pdb_id, chains_1, chains_2 = pdb.split("_") 

    chain2seq = dict()
    for c in chains_1:
        chain2seq[c] = deepcopy(all_sequences_wt[pdb_id][c])
    for c in chains_2:
        chain2seq[c] = deepcopy(all_sequences_wt[pdb_id][c])
    
    ### mutated types
    mutations = row_dict["Mutation(s)_cleaned"]
    mutations = mutations.split(",")
    chain2mutations = defaultdict(list)
    for mut in mutations:
        wt = mut[0]
        chain_id = mut[1]
        pos = int(mut[2:-1])
        mt = mut[-1]
        seq = chain2seq[chain_id]
        # assert seq[pos-1] == wt  
        seq[pos-1] = mt
        chain2seq[chain_id] = seq
        chain2mutations[chain_id].append(f"{wt}{pos}{mt}")

    interactor_1 = [] # could have multiple chains
    interactor_2 = []
    for chain_id in chains_1:
        seq = "".join(chain2seq[chain_id])
        id = f"{pdb_id}_{chain_id}"
        if len(chain2mutations[chain_id]) > 0:
            id = f"{id}_mut_{'_'.join(chain2mutations[chain_id])}"
        interactor_1.append({
            "pdb_chain": f"{pdb_id}_{chain_id}",
            "mutations": chain2mutations[chain_id],
            "seq": seq,
            # "id": id
        })
        seq2ids[seq].append(id)
    
    for chain_id in chains_2:
        seq = "".join(chain2seq[chain_id])
        id = f"{pdb_id}_{chain_id}"
        if len(chain2mutations[chain_id]) > 0:
            id = f"{id}_mut_{'_'.join(chain2mutations[chain_id])}"
        interactor_2.append({
            "pdb_chain": f"{pdb_id}_{chain_id}",
            "mutations": chain2mutations[chain_id],
            "seq": seq,
            # "id": id
        }) 
        seq2ids[seq].append(id)
    
    interactor_1 = sorted(interactor_1, key=lambda x: x["seq"])
    interactor_2 = sorted(interactor_2, key=lambda x: x["seq"])

    if tuple([x["seq"] for x in interactor_1]) > tuple([x["seq"] for x in interactor_2]):
        interactor_1, interactor_2 = interactor_2, interactor_1
        # row_dict["Protein 1"], row_dict["Protein 2"] = row_dict["Protein 2"], row_dict["Protein 1"]

    # print((tuple(x["seq"] for x in interactor_1), tuple(x["seq"] for x in interactor_2)))
    # break
    pdb2records[pdb][(tuple(x["seq"] for x in interactor_1), tuple(x["seq"] for x in interactor_2))].append(
        {
            "interactors": [interactor_1, interactor_2], 
            "interaction_attr": row_dict,
            "wt": False
        })

# print(len(pdb2records)) # 
# for pdb in pdb2records:
#     for seq_pairs in pdb2records[pdb]:
#         if len(pdb2records[pdb][seq_pairs]) > 1:
#             print(f"#pdb: {pdb}")
#             print(pdb2records[pdb][seq_pairs])
#             exit()

seq2unique_id = dict()
with open(output_dir / "unique_sequences.fasta", "w") as fout:
    for seq in seq2ids:
        all_ids = list(set(seq2ids[seq]))
        all_ids.sort()
        seq2unique_id[seq] = all_ids[0]
        # fout.write(f">{all_ids[0]} {' '.join(all_ids[1:])}\n{seq}\n") # no need to save this?!
        fout.write(f">{all_ids[0]}\n{seq}\n")

unf_entries = 0
conflict_binding = 0
records = []
for pdb in pdb2records:
    for seqs in pdb2records[pdb]:
        interactors_1, interactors_2 = seqs
        interactors = []
        for seq_1 in interactors_1:
            interactors.append({"seq": seq_1, "id": seq2unique_id[seq_1], "role": 0})
        for seq_2 in interactors_2:
            interactors.append({"seq": seq_2, "id": seq2unique_id[seq_2], "role": 1})
            
        experiments = pdb2records[pdb][seqs]

        attributes = []    
        affinities = [] # just for check!
        for experiment in experiments:
            interactor_attrs = []
            for interactor_idx, interactor in enumerate(experiment["interactors"]):
                for chain in interactor:
                    chain = deepcopy(chain)
                    seq = chain.pop("seq")
                    chain["interactor_idx"] = interactor_idx
                    interactor_attrs.append(chain)
            attributes.append(
                {
                    "interaction_attr": experiment["interaction_attr"], 
                    "interactor_attr": interactor_attrs,
                    "wt": experiment["wt"]
                })        
            if not experiment["wt"]:
                affinities.append((experiment["interaction_attr"]["Method"], 
                                experiment["interaction_attr"]["Affinity_mut_parsed"],
                                experiment["interaction_attr"]["Temperature"]))
            else:
                affinities.append((experiment["interaction_attr"]["Method"], 
                                experiment["interaction_attr"]["Affinity_wt_parsed"],
                                experiment["interaction_attr"]["Temperature"]))
        
        unf = False
        bindings = []
        for attribute in attributes:
            interaction_attr = attribute["interaction_attr"]
            if attribute["wt"]:
                affinity_parsed = interaction_attr["Affinity_wt_parsed"]
                affinity = interaction_attr["Affinity_wt (M)"]
            else:
                affinity_parsed = interaction_attr["Affinity_mut_parsed"]
                affinity = interaction_attr["Affinity_mut (M)"]
            
            binding = None
            if affinity == "unf":
                unf = True
            elif affinity == "n.b" or affinity == "n.b.":
                binding = False
            elif affinity_parsed >= 5e-6:
                binding = False
            else:
                binding = True
            # elif affinity_parsed >= 1e-4:
            #     binding = 0
            # elif affinity_parsed >= 1e-5:
            #     binding = 1
            # else:
            #     binding = 2
            if binding is not None:
                bindings.append(binding)

        if unf: # skip unf
            unf_entries += 1
            continue
        bindings_counter = Counter(bindings).most_common()

        # if len(bindings_counter) > 1:
            # print(f"warning: conflict binding affintiy: {bindings_counter} from {affinities}")
                    
        if bindings_counter[0][1] / len(bindings) > 0.5:
            binding = bindings_counter[0][0]
        else:
            print(f"warning: conflict binding affintiy: {bindings_counter} from {affinities}")
            conflict_binding += 1
            continue
            # if 0 in bindings and 2 in bindings:
            #     conflict_binding += 1
            #     continue
            # bindings.sort()
            # binding = bindings[0]

        records.append({
            "interactors": interactors, 
            "attributes": attributes, 
            "assay_id": f"skempi_{pdb}",
            "binding": binding
            })
        
        # "Affinity_mut (M)": "n.b", "Affinity_mut_parsed": ""

print(f"unf_entries: {unf_entries}")
print(f"conflict_binding: {conflict_binding}")
print(f"records: {len(records)}")
json.dump(records, open(output_dir / "skempi_v2_by_pdb.json", "w"))