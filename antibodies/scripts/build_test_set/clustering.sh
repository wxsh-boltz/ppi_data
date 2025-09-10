

root_dir="/data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/sabdab_summary_all_1vh1vl.test_vh_0.7_vl_0.7_an_0.7"
output_dir="$root_dir/clusterRes/"
mkdir -p $output_dir

mmseqs easy-cluster $root_dir/antigen_seqs.fasta $output_dir/antigen tmp --min-seq-id 0.7 -c 0.25 --cov-mode 1
mmseqs easy-cluster $root_dir/vh_seqs.fasta $output_dir/vh tmp --min-seq-id 0.7 -c 0.25 --cov-mode 1
mmseqs easy-cluster $root_dir/vl_seqs.fasta $output_dir/vl tmp --min-seq-id 0.7 -c 0.25 --cov-mode 1

python 0_calc_edit_distance_for_peptides.py \
    --query_fasta_path $root_dir/antigen_seqs.fasta \
    --target_fasta_path $root_dir/antigen_seqs.fasta \
    --saving_path_alignment $root_dir/antigen_peptide_TO_antigen_peptide.m8 \
    --saving_path_peptides $root_dir/antigen_peptide.txt