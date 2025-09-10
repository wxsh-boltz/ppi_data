cov_mode=0

e=0.1 # 10000000

ANTIGEN_PATH="/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/antigen_seqs.fasta"
RCSB_SEQ_PATH="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/pdb_05_08_2025_wx/rcsb_boltz2_train_seq_boltz_extended_and_gemmi.fasta"
RCSB_VH_PATH="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/pdb_05_08_2025_wx/anarci/vh.fasta"
RCSB_VL_PATH="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/pdb_05_08_2025_wx/anarci/vl.fasta"

VH_PATH="/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/vh_seqs.fasta"
VL_PATH="/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/vl_seqs.fasta"

alignment_output_dir="/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/alnRes"
mkdir -p $alignment_output_dir

# antigen TO rcsb
cov_mode="1"
OUTPUR_PATH="$alignment_output_dir/antigen_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8"
mmseqs easy-search $ANTIGEN_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e

cov_mode="2"
OUTPUR_PATH="$alignment_output_dir/antigen_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8"
mmseqs easy-search $ANTIGEN_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e

# VH TO rcsb (vh)
# cov_mode="1"
# OUTPUR_PATH="$alignment_output_dir/vh_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8"
# mmseqs easy-search $VH_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 -e 10000000 --spaced-kmer-mode 0
OUTPUR_PATH="$alignment_output_dir/vh_seqs_TO_rcsb_boltz2_train_vhs_cov_mode_0.m8"
mmseqs easy-search $VH_PATH $RCSB_VH_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e

# VL to rcsb
# cov_mode="1"
# OUTPUR_PATH="$alignment_output_dir/vl_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8"
# mmseqs easy-search $VL_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 -e 10000000 --spaced-kmer-mode 0
OUTPUR_PATH="$alignment_output_dir/vl_seqs_TO_rcsb_boltz2_train_vls_cov_mode_0.m8"
mmseqs easy-search $VL_PATH $RCSB_VL_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e


#### calculate this just to figure which sequences are not processed by the mmseqs
# target TO target
OUTPUR_PATH="$alignment_output_dir/antigen_seqs_TO_antigen_seqs.m8"
mmseqs easy-search $ANTIGEN_PATH $ANTIGEN_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 1.0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e

# VH TO VH
OUTPUR_PATH="$alignment_output_dir/vh_seqs_TO_vh_seqs.m8"
mmseqs easy-search $VH_PATH $VH_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 1.0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e

# VL TO VL
OUTPUR_PATH="$alignment_output_dir/vl_seqs_TO_vl_seqs.m8"
mmseqs easy-search $VL_PATH $VL_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 1.0 --max-seqs 10000000 --spaced-kmer-mode 0 -e $e