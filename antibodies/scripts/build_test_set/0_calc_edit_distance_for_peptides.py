from Bio import SeqIO
from collections import defaultdict
import Levenshtein
from tqdm import tqdm
import re, sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Peptide extraction and alignment configuration")

    parser.add_argument("--max_peptide_length", type=int, default=14,
                        help="Maximum length of peptides to extract")
    parser.add_argument("--max_edit_distance", type=int, default=3,
                        help="Maximum edit distance allowed for peptide alignment")

    parser.add_argument("--target_fasta_path", type=str, # required=True, 
                        default="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/pdb_05_08_2025_wx/rcsb_boltz2_train_seq_boltz_extended_and_gemmi.fasta",
                        help="Path to the target FASTA file")

    parser.add_argument("--query_fasta_path", type=str, required=True,
                        help="Path to the query FASTA file")

    parser.add_argument("--saving_path_alignment", type=str, required=True,
                        help="Path to save the alignment results")

    parser.add_argument("--saving_path_peptides", type=str, required=True,
                        help="Path to save the peptide results")

    return parser.parse_args()


def custom_edit_distance(s1, s2, wildcard='X'):
    len1, len2 = len(s1), len(s2)
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]

    for i in range(len1 + 1):
        dp[i][0] = i
    for j in range(len2 + 1):
        dp[0][j] = j

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            c1, c2 = s1[i - 1], s2[j - 1]
            if c1 == c2 or c1 == wildcard or c2 == wildcard:
                cost = 0
            else:
                cost = 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,      # deletion
                dp[i][j - 1] + 1,      # insertion
                dp[i - 1][j - 1] + cost  # substitution
            )

    return dp[len1][len2]

### peptides: <= 14 AAs
# max_peptide_length = 14
# max_edit_distance = 3

# target_fasta_path = "/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/pdb_05_08_2025_wx/rcsb_boltz2_train_seq_boltz_extended_and_gemmi.fasta"

# # snac:
# # query_fasta_path = "/data/rbg/users/wxsh/esm3/data/antibodies/SNAC/processed/snac_antigen_seqs.fasta"
# # saving_path_alignment = "/data/rbg/users/wxsh/esm3/data/antibodies/SNAC/processed/alnRes/antigen_peptide_TO_rcsb_boltz2_train_seqs.m8"
# # saving_path_peptides = "/data/rbg/users/wxsh/esm3/data/antibodies/SNAC/processed/antigen_peptide.txt"

# # SabDab
# query_fasta_path = "/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/antigen_seqs.fasta"
# saving_path_alignment = "/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/alnRes/antigen_peptide_TO_rcsb_boltz2_train_seqs.m8"
# saving_path_peptides = "/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/antigen_peptide.txt"

args = parse_args()

query_peptides = defaultdict(list)
for record in SeqIO.parse(args.query_fasta_path, "fasta"):
    if len(record.seq) > args.max_peptide_length:
        continue
    # query_peptides.append((record.id, str(record.seq)))
    # query_seq = re.sub(r'[a-z]', 'X', str(record.seq))
    query_seq = str(record.seq).upper()
    query_peptides[query_seq].append(record.id)
print(f"Query peptides: {len(query_peptides)}")

with open(args.saving_path_peptides, "w") as fout:
    for seq, ids in query_peptides.items():
        # assert len(ids) == 1, f"{seq}, {ids}"
        for id in ids:
            fout.write(f"{id}\n")

target_peptides = defaultdict(list)
for record in SeqIO.parse(args.target_fasta_path, "fasta"):
    if len(record.seq) > args.max_peptide_length + args.max_edit_distance: # at most insert max_edit_distance AAs
        continue
    # target_seq = re.sub(r'[a-z]', 'X', str(record.seq))
    target_seq = str(record.seq).upper()
    # target_seq = target_seq.upper()
    target_peptides[target_seq].append(record.id)
    # target_peptides.append((record.id, str(record.seq)))
print(f"Target peptides: {len(target_peptides)}")

aligments = []
for query_seq in tqdm(query_peptides):
    for target_seq in target_peptides:
        distance = Levenshtein.distance(query_seq, target_seq)
        # distance = custom_edit_distance(query_seq, target_seq)
        # if "X" in target_seq and distance == 4:
        #     print(target_seq)
        #     print(query_seq)
        #     print(custom_edit_distance(query_seq, target_seq))
        #     print(Levenshtein.distance(query_seq, target_seq))
        #     # exit()
        # elif "X" not in target_seq and "X" not in query_seq:
        #     assert custom_edit_distance(query_seq, target_seq) == distance
        if distance <= args.max_edit_distance:
            for query_id in query_peptides[query_seq]:
                for target_id in target_peptides[target_seq]:
                    aligments.append((query_id, target_id, distance))


print("aligments", len(aligments))
with open(args.saving_path_alignment, "w") as fout:
    for query_id, target_id, distance in aligments:
        fout.write(f"{query_id}\t{target_id}\t{distance}\n")