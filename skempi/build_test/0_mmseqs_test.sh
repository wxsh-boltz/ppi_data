RCSB_SEQ_PATH="/data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/pdb_05_08_2025_wx/rcsb_boltz2_train_seq_boltz_extended_and_gemmi.fasta"
SEQ_PATH="/data/rbg/users/wxsh/esm3/data/skempi/processed/unique_sequences.fasta"

data_root_dir="/data/rbg/users/wxsh/esm3/data/skempi/processed"
alignment_output_dir="$data_root_dir/alnRes"
mkdir -p $alignment_output_dir

cov_mode="1"
OUTPUR_PATH="$alignment_output_dir/unique_sequences_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8"
mmseqs easy-search $SEQ_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 -e 0.1 --spaced-kmer-mode 0 

cov_mode="2"
OUTPUR_PATH="$alignment_output_dir/unique_sequences_TO_rcsb_boltz2_train_seqs_cov_mode_$cov_mode.m8"
mmseqs easy-search $SEQ_PATH $RCSB_SEQ_PATH $OUTPUR_PATH tmp --cov-mode $cov_mode -c 0.25 --min-seq-id 0 --max-seqs 10000000 -e 0.1 --spaced-kmer-mode 0

OUTPUR_PATH="$alignment_output_dir/unique_sequences_TO_unique_sequences.m8"
mmseqs easy-search $SEQ_PATH $SEQ_PATH $OUTPUR_PATH tmp --cov-mode 0 -c 0.25 --min-seq-id 1.0 --max-seqs 10000000 -e 0.1 --spaced-kmer-mode 0 

python /data/rbg/users/wxsh/esm3/data/antibodies/scripts/build_test_set/0_calc_edit_distance_for_peptides.py \
    --max_edit_distance 0 \
    --query_fasta_path $SEQ_PATH \
    --saving_path_alignment "$alignment_output_dir/peptide_sequences_TO_rcsb_boltz2_train_seqs.m8" \
    --saving_path_peptides "$alignment_output_dir/peptide_sequences.fasta"
