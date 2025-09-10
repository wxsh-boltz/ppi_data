
Code for building antibody test / validation set and excluding the leakage from RCSB training set.

## Step 1: calculate sequence similarity

First, we need to calculate the sequence identity with RCSB train using `mmseqs easy-search`.
Because mmseqs cannot handle short sequences, we will also calculate the edit distance for short antigens (VH/VL usually won't be very short).

```
bash build_test_set/0_mmseqs_search_sabdab_new.sh

RCSB_SEQ_FASTA_PATH="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/pdb_05_08_2025_wx/rcsb_boltz2_train_seq_boltz_extended_and_gemmi.fasta"
ANTIGEN_SEQ_FASTA_PATH="/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/antigen_seqs.fasta"
SAVING_PEPTIDE_PATH="/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/antigen_peptide.txt"
SAVING_ANTIGEN_PEPTIDE_ALN="/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/alnRes/antigen_peptide_TO_rcsb_boltz2_train_seqs.m8"

python build_test_set/0_calc_edit_distance_for_peptides.py --target_fasta_path $RCSB_SEQ_FASTA_PATH --query_fasta_path $ANTIGEN_SEQ_FASTA_PATH --saving_path_alignment $SAVING_ANTIGEN_PEPTIDE_ALN --saving_path_peptides $SAVING_PEPTIDE_PATH
```

## Step 2: exclude leakage

```bash

root_dir="/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed"

python 1_split_test_antibody.py \
    --input_path $root_dir/sabdab_summary_all_1vh1vl.json \
    --vh_aln_paths $root_dir/alnRes/vh_seqs_TO_rcsb_boltz2_train_vhs_cov_mode_0.m8 \
    --vl_aln_paths $root_dir/alnRes/vl_seqs_TO_rcsb_boltz2_train_vls_cov_mode_0.m8 \
    --antigen_aln_paths $root_dir/alnRes/antigen_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_1.m8  $root_dir/alnRes/antigen_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_2.m8 \
    --vh_self_aln_path $root_dir/alnRes/vh_seqs_TO_vh_seqs.m8 \
    --vl_self_aln_path $root_dir/alnRes/vl_seqs_TO_vl_seqs.m8 \
    --antigen_self_aln_path $root_dir/alnRes/antigen_seqs_TO_antigen_seqs.m8 \
    --vh_sequence_identity_threshold 0.7 \
    --vl_sequence_identity_threshold 0.7 \
    --antigen_sequence_identity_threshold 0.7 \
    --peptide_antigen_aln_path $root_dir/alnRes/antigen_peptide_TO_rcsb_boltz2_train_seqs.m8 \
    --peptide_antigen_list_path $root_dir/antigen_peptide.txt \
    --assay_id "sabdab_test"
```

## Step 3: build decoys

Hard code for sabdab, but you can easily change to other datasets.

```bash
# 1. clustering vh/vl/antigens
bash clustering.sh

# 2. build decoys
python 3_build_decoy_new.py

```

## Step 4: subsample

The following scripts are used to subsample the test set. 

`4_downsampling_random.py`: randomly sample $N$ interactions for each antigen/binding label. 

`4_downsampling_random_neg.py`: keep all the positives, and subsample the negatives to $N$.
