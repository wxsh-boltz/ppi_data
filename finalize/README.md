For constructing the training set and excluding the leakage from ALL the testing/validation set.

### Step 1: extract new targets from test/val sets

First, I exctract the target (or antigen) from antibody or peptide test/val sets that do not occur in any complexes in the RCSB training set.

Run `python 0_collect_new_targets.py`

### Step 2: exclude leakage from training set

Exclude the leakage for all training set from all the testing/validation sets. Config file: `datasets.yaml`.

Run `python 1_split_training_set.py`