source activate anarci
conda activate anarci

cd /data/rbg/users/wxsh/ANARCI

echo "Run ANARCI for hchain"
ANARCI -i /data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/hchain_seqs.fasta -o /data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/hchain_seqs.anarci --csv

echo "Run ANARCI for lchain"
ANARCI -i /data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/lchain_seqs.fasta -o /data/rbg/users/wxsh/esm3/data/antibodies/sabdab/processed/lchain_seqs.anarci --csv


# /data/rbg/users/wxsh/esm3/data/boltz2_training/rcsb/polymers_after_20230601/rcsb_2prot_after_20230601.test_si_0.7/seqs.fasta