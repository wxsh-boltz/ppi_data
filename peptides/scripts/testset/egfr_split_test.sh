root_dir="/data/rbg/users/wxsh/esm3/data/peptides/egfr_competition_2/"
python 1_split_testset.py \
    --target_align_path_mode_2 "$root_dir/alnRes/target_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_2.m8" \
    --target_align_path_mode_1 "$root_dir/alnRes/target_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_1.m8" \
    --binder_align_path "$root_dir/alnRes/binder_seqs_TO_rcsb_boltz2_train_seqs_cov_mode_1.m8" \
    --binder_self_align_path "$root_dir/alnRes/binder_seqs_TO_binder_seqs.m8" \
    --target_self_align_path "$root_dir/alnRes/target_seqs_TO_target_seqs.m8" \
    --data_path $root_dir/egfr_competition_2.json \
    --binder_si_threshold 0.7 --target_si_threshold 0.7