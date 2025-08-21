# data_name="mini_proteins_chevalier_nature_2017"
# TARGET_PATH="/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/"$data_name"_target_seqs.fasta"
# BINDER_PATH="/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/"$data_name"_binder_seqs.fasta"
# RCSB_SEQ_PATH="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/rcsb_boltz2_train_seqs.fasta"
# alignment_output_dir="/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/alnRes"
# mkdir $alignment_output_dir

# # target TO rcsb
# cov_mode=1
# OUTPUR_PATH="$alignment_output_dir/target_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8"
# mmseqs easy-search $TARGET_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 

# cov_mode=2
# OUTPUR_PATH="$alignment_output_dir/target_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8"
# mmseqs easy-search $TARGET_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 

# # binders TO rscb
# OUTPUR_PATH="$alignment_output_dir/binder_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_1.m8"
# mmseqs easy-search $BINDER_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode 1 -c 0.25 --min-seq-id 0 --max-seqs 10000000 

# # target TO target
# OUTPUR_PATH="$alignment_output_dir/target_seqs_TO_target_seqs.m8"
# mmseqs easy-search $TARGET_PATH $TARGET_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 1.0 --max-seqs 10000000

# # binder TO binder
# OUTPUR_PATH="$alignment_output_dir/binder_seqs_TO_binder_seqs.m8"
# mmseqs easy-search $BINDER_PATH $BINDER_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 1.0 --max-seqs 10000000 

e=0.1
RCSB_SEQ_PATH="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/pdb_05_08_2025_wx/rcsb_boltz2_train_seq_boltz_extended_and_gemmi.fasta"

### 2022
TARGET_PATH="/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/mini_proteins_cao_nature_2022_target_seqs.fasta"
BINDER_PATH="/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/mini_proteins_cao_nature_2022_binder_seqs.fasta"
alignment_output_dir="/data/rbg/users/wxsh/esm3/data/peptides/baker_2022/processed_merged/alnRes"

# ### 2017
TARGET_PATH="/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/mini_proteins_chevalier_nature_2017_target_seqs.fasta"
BINDER_PATH="/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/mini_proteins_chevalier_nature_2017_binder_seqs.fasta"
alignment_output_dir="/data/rbg/users/wxsh/esm3/data/peptides/baker_2017/processed_merged/alnRes"


mkdir -p $alignment_output_dir

# target TO rcsb
cov_mode="1"
mmseqs easy-search $TARGET_PATH $RCSB_SEQ_PATH "$alignment_output_dir/target_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8" tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e 

cov_mode="2"
mmseqs easy-search $TARGET_PATH $RCSB_SEQ_PATH "$alignment_output_dir/target_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8" tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e 

# binder TO rcsb
cov_mode="1"
OUTPUR_PATH="$alignment_output_dir/binder_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8"
mmseqs easy-search $BINDER_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e 

# target TO target
OUTPUR_PATH="$alignment_output_dir/target_seqs_TO_target_seqs.m8"
mmseqs easy-search $TARGET_PATH $TARGET_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 1.0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e 

# binder TO binder
OUTPUR_PATH="$alignment_output_dir/binder_seqs_TO_binder_seqs.m8"
mmseqs easy-search $BINDER_PATH $BINDER_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 1.0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e 