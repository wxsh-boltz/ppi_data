import pandas as pd
import json
from collections import defaultdict
import re

direct_interaction_id = "MI:0407"
input_path = "../BIOGRID-MV-Physical-4.4.245.mitab.txt"
# output_path = input_path.replace("_all_interactions.csv", "_direct_interactions.csv")
output_path = input_path.replace(".mitab.txt", ".direct_interactions.csv")
# interaction_type_column = "intact_types"
interaction_type_column = "Interaction Types"

def read_intact_ontology(path="../../uniprot/psi-mi.csv"):
    # def _get_
    node2parent_node = dict()
    node2babies = defaultdict(list)
    df = pd.read_csv(path)
    for node, parent_node in zip(df["id"], df["subclass_of"]):
        node2parent_node[node] = parent_node
        node2babies[parent_node].append(node)

    node2offspring_nodes = defaultdict(list)
    def _get_offsprings(node, l):
        for baby in node2babies.get(node, []):
            l.append(baby)
            _get_offsprings(baby, l)
        return l

    for node in node2babies:
        node2offspring_nodes[node] = _get_offsprings(node, [])
    return node2offspring_nodes

node2offspring_nodes = read_intact_ontology()
direct_interactions = node2offspring_nodes["MI:0407"] + ["MI:0407"]
direct_interactions = set(direct_interactions)
print(f"direct_interactions classes: {len(direct_interactions)}")

def extract_mi_id(s):
    match = re.search(r'MI:\d{4}', s)
    if match:
        extracted_id = match.group()
        return extracted_id
    else:
        print("No MI:xxxx ID found")

if input_path.endswith(".mitab.txt"):
    df = pd.read_csv(input_path, sep="\t")
elif input_path.endswith(".csv"):
    df = pd.read_csv(input_path)
else:
    raise ValueError(f"Not support reading {input_path}")

da_mask = []
for intact_type in df[interaction_type_column]:
    da = False
    int_id = extract_mi_id(intact_type)
    if int_id in direct_interactions:
        da = True
    da_mask.append(da)
df_direct = df[da_mask]
print(f"Input {len(df)}, direct interactions: {len(df_direct)}")
df_direct.to_csv(output_path, index=False)


