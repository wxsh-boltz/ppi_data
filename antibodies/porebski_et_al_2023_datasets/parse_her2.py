import pandas as pd
import numpy as np
from collections import defaultdict
import json
from abnumber import Chain
from pathlib import Path
from tqdm import tqdm

# Long
HER2_seq="MELAALCRWGLLLALLPPGAASTQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPGGLRELQLRSLTEILKGGVLIQRNPQLCYQDTILWKDIFHKNNQLALTLIDTNRSRACHPCSPMCKGSRCWGESSEDCQSLTRTVCAGGCARCKGPLPTDCCHEQCAAGCTGPKHSDCLACLHFNHSGICELHCPALVTYNTDTFESMPNPEGRYTFGASCVTACPYNYLSTDVGSCTLVCPLHNQEVTAEDGTQRCEKCSKPCARVCYGLGMEHLREVRAVTSANIQEFAGCKKIFGSLAFLPESFDGDPASNTAPLQPEQLQVFETLEEITGYLYISAWPDSLPDLSVFQNLQVIRGRILHNGAYSLTLQGLGISWLGLRSLRELGSGLALIHHNTHLCFVHTVPWDQLFRNPHQALLHTANRPEDECVGEGLACHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPINCTHSCVDLDDKGCPAEQRASPLTSIISAVVGILLVVVLGVVFGILIKRRQQKIRKYTMRRLLQETELVEPLTPSGAMPNQAQMRILKETELRKVKVLGSGAFGTVYKGIWIPDGENVKIPVAIKVLRENTSPKANKEILDEAYVMAGVGSPYVSRLLGICLTSTVQLVTQLMPYGCLLDHVRENRGRLGSQDLLNWCMQIAKGMSYLEDVRLVHRDLAARNVLVKSPNHVKITDFGLARLLDIDETEYHADGGKVPIKWMALESILRRRFTHQSDVWSYGVTVWELMTFGAKPYDGIPAREIPDLLEKGERLPQPPICTIDVYMIMVKCWMIDSECRPRFRELVSEFSRMARDPQRFVVIQNEDLGPASPLDSTFYRSLLEDDDMGDLVDAEEYLVPQQGFFCPDPAPGAGGMVHHRHRSSSTRSGGGDLTLGLEPSEEEAPRSPLAPSEGAGSDVFDGDLGMGAAKGLQSLPTHDPSPLQRYSEDPTVPLPSETDGYVAPLTCSPQPEYVNQPDVRPQPPSPREGPLPAARPAGATLERPKTLSPGKNGVVKDVFAFGGAVENPEYLTPQGGAAPQPHPPPAFSPAFDNLYYWDQDPPERGAPPSTFKGTPTAENPEYLGLDVPV"
HER2_id="P04626"

# HER2_pdb = "TQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPGGLRELQLRSLTEILKGGVLIQRNPQLCYQDTILWKDIFHKNNQLALTLIDTNRSRACHPCSPMCKGSRCWGESSEDCQSLTRTVCAGGCARCKGPLPTDCCHEQCAAGCTGPKHSDCLACLHFNHSGICELHCPALVTYNTDTFESMPNPEGRYTFGASCVTACPYNYLSTDVGSCTLVCPLHNQEVTAEDGTQRCEKCSKPCARVCYGLGMEHLREVRAVTSANIQEFAGCKKIFGSLAFLPESFDGDPASNTAPLQPEQLQVFETLEEITGYLYISAWPDSLPDLSVFQNLQVIRGRILHNGAYSLTLQGLGISWLGLRSLRELGSGLALIHHNTHLCFVHTVPWDQLFRNPHQALLHTANRPEDECVGEGLACHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPIN"
# len(HER2_uniprot), len(HER2_pdb)

dataset_name2path = {
    "porebski_her2_affmat": "3.HER2_scFv/fc079_fc080_her2_affmat_paired_translated.csv",
    "porebski_her2_ml": "3.HER2_scFv/fc081_fc082_her2_ml_lib_paired_translated.csv"
}
dataset_name = "porebski_her2_affmat"
path = dataset_name2path[dataset_name]

saving_path = Path("processed") / dataset_name
saving_path.mkdir(exist_ok=True, parents=True)

# matched with Talip's data
vl = "QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYGHTNRPAGVPDRFSGSKSGTSASLAISGFRSEDEADYYCASWDYTLSGWVFGGGTKLTVL"
# vl_check = str(Chain(vl, scheme='imgt'))
# print(vl_check == vl)

df=pd.read_csv(path)

unique_seq_scores = defaultdict(list)
for _, row in tqdm(df.iterrows(), total=len(df)):
    row_list = row.tolist()
    seq = row_list[1]
    targets = row_list[2:]
    
    if "PPPPPPP" in seq:
        continue
    if "KKKKKKK" in seq:
        continue
    if "GGGGGGG" in seq:
        continue
    if "*" in seq: # shoud we?
        continue
        
    if targets[0] >= 0.0:
        unique_seq_scores[seq].append(targets)
    
dataset = []

label_counting = defaultdict(int)
for i, seq in tqdm(enumerate(unique_seq_scores), total=len(unique_seq_scores)):
    targets = np.array(unique_seq_scores[seq])
    conc = np.max(targets, axis=0)
    tgt_condition = conc[8]

    # Each sequence was given a label based on its FImean in the 5 minute wash condition
    if tgt_condition < 150.0: ## -> parent seqs
        peak_tgt = 0
    if tgt_condition >= 150.0 and tgt_condition < 250.0:
        peak_tgt = 1
    if tgt_condition >= 250.0: ## 3 sigma 
        peak_tgt = 2

    label_counting[peak_tgt] += 1

    # Full VH sequence:
    seq_composed = "QVQLVQSGAEVKKPGESLKISCKGSGYSFTSYWIAWVRQMPGKGLEYMGLIYPGDSDTKYSPSFQGQVTISVDKSVSTAYLQWSSLKPSDSAVYFCAR%sGQGTLVTVSS" % (seq)
    dataset.append([seq_composed, peak_tgt])

    # vh_check = str(Chain(seq_composed, scheme='imgt'))
    # assert vh_check == seq_composed
    # print(vh_check == seq_composed)
    # exit()

print(label_counting)

all_vh_seqs = [x[0] for x in dataset]
all_vh_seqs.sort()
vh_seq2id = {}
for i, seq in enumerate(all_vh_seqs):
    vh_seq2id[seq] = f"{dataset_name}_vh_{i}"

print(f"Number of VHs", len(vh_seq2id))
with open(saving_path / f"vh_seqs.fasta", "w") as fout:
    for seq in vh_seq2id:
        fout.write(f">{vh_seq2id[seq]}\n{seq}\n")

with open(saving_path / f"vl_seqs.fasta", "w") as fout:
    fout.write(f">anti_her2_g98a_vl\n{vl}\n")

with open(saving_path / f"antigen_seqs.fasta", "w") as fout:
    fout.write(f">{HER2_id}\n{HER2_seq}\n")

records = []
for seq, label in dataset:
    interactors = []
    interactors.append({"seq": seq, "id": vh_seq2id[seq], "role": "vh"})
    interactors.append({"seq": vl, "id": "anti_her2_g98a_vl", "role": "vl"})
    interactors.append({"seq": HER2_seq, "id": HER2_id, "role": "antigen"})
    records.append({
        "interactors": interactors, 
        "binding": label, 
        "assay_id": dataset_name})
print("records:", len(records))
json.dump(records, open(saving_path / f"{dataset_name}.json", "w"))