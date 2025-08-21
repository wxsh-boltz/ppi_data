source activate anarci
conda activate anarci

cd /data/rbg/users/wxsh/ANARCI

echo "Run ANARCI for hchain"
ANARCI -i /data/rbg/users/wxsh/esm3/data/antibodies/AVIDa-hIL6/AVIDa-hIL6_vhh_seqs.fasta -o /data/rbg/users/wxsh/esm3/data/antibodies/AVIDa-hIL6/AVIDa-hIL6_vhh_seqs.anarci.csv --csv