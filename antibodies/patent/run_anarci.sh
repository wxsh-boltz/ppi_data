
source activate anarci

cd /data/rbg/users/wxsh/ANARCI

echo "Run ANARCI for hchain"
ANARCI -i /data/rbg/users/wxsh/esm3/data/antibodies/patent/heavy_sequences.fasta -o /data/rbg/users/wxsh/esm3/data/antibodies/patent/heavy_sequences.anarci --csv

echo "Run ANARCI for lchain"
ANARCI -i /data/rbg/users/wxsh/esm3/data/antibodies/patent/light_sequences.fasta -o /data/rbg/users/wxsh/esm3/data/antibodies/patent/light_sequences.anarci --csv