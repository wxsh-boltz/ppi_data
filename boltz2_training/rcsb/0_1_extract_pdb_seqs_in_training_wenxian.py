import json
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict

all_seqs_path = "pdb_05_08_2025_wx/sequences.json"
date_cutoff = np.datetime64("2023-06-01")
saving_path_prefix = "pdb_05_08_2025_wx/rcsb_boltz2_train"

# {'id': '2afr',
#   'released': '2006-04-04',
#   'proteins': [['2afr',
#     0,
#     '1',
#     'A1',
#     'MNDMRQITNLGRNIENKSFSIIDEEAGPHSFAQEEWEVVRRIIHATADFDYKNITKIHPQAIDSGIQALKKGCPIVCDVQMILSGLNPERLKVYGCKTYCFISDEDVIENAKRKNSTRAIESIQKANSFNLLNESIIVIGNAPTALLEIEKLIRQEGIKPALIVGVPVGFVSAKESKESILKLEYYNVTSIPYILTMGRKGGSTIAVAILHALLLLSSKRGERLEHHHHHH']],
#   'rnas': [],
#   'dnas': [],
#   'ligands': []})

dataset = json.load(open(all_seqs_path))
print(f"Read {len(dataset)} data")
id2seq = dict()

for entry in dataset:
    if np.datetime64(entry["released"]) > date_cutoff: # skip 
        continue
    for protein in entry["proteins"]:
        pdb_id, chain_id, _, _, seq = protein
        assert f"{pdb_id}_{chain_id}" not in id2seq, protein

        for seq_from in seq:
            if seq_from not in id2seq:
                id2seq[seq_from] = dict()
            id2seq[seq_from][f"{pdb_id}_{chain_id}"] = seq[seq_from]

for seq_from in id2seq:
    print(f"Read {len(id2seq)} seqs ({len(id2seq[seq_from])})")

for seq_from in id2seq:
    with open(f"{saving_path_prefix}_{seq_from}.fasta", "w") as fout:
        for seq_id in id2seq[seq_from]:
            fout.write(f">{seq_id}\n{id2seq[seq_from][seq_id]}\n")