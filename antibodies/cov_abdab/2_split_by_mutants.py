import json
from collections import defaultdict, Counter
from pathlib import Path
input_data = "/data/rbg/users/wxsh/esm3/data/antibodies/cov_abdab/processed/CoV-AbDab_080224.json"
output_dir = Path("/data/rbg/users/wxsh/esm3/data/antibodies/cov_abdab/processed/CoV-AbDab_080224_split_by_variants")
output_dir.mkdir(exist_ok=True, parents=True)

all_data = json.load(open(input_data))
print(f"All: {len(all_data)}, {Counter([d['binding'] for d in all_data])}")

# print(all_data[0])

variant2data = defaultdict(list)
training_variants = [
    "SARS-CoV2_WT",
    "SARS-CoV1",
    "SARS-CoV2_Alpha",
    "SARS-CoV2_Beta",
    # "SARS-CoV2_Gamma",
    # "SARS-CoV2_Delta",
    ]

for data in all_data:
    antigen = [interactor["id"] for interactor in data["interactors"] if interactor["role"] == "antigen"][0]
    variant = antigen.split("_")[:2]
    if variant[0] == "SARS-CoV1":
        variant = "SARS-CoV1"
    else:
        variant = variant[0] + "_" + variant[1].split("-")[0]
    variant2data[variant].append(data)

    
print(len(variant2data))

training = []
testing = []
testing_variants = set()
for v in variant2data:
    print(v, len(variant2data[v]))
    if v in training_variants:
        training.extend(variant2data[v])
    else:
        testing_variants.add(v)
        testing.extend(variant2data[v])
print(f"Train: {len(training)}, {Counter([d['binding'] for d in training])}")
print(f"Test: {len(testing)}, {Counter([d['binding'] for d in testing])}")
print(testing_variants)
with open(output_dir / "variants.csv", "w") as fout:
    fout.write(f"Variant,Split\n")
    for v in testing_variants:
        fout.write(f"{v},Test\n")
    for v in training_variants:
        fout.write(f"{v},Train\n")

json.dump(training, open(output_dir / "train.json", "w"))
json.dump(testing, open(output_dir / "test.json", "w"))
