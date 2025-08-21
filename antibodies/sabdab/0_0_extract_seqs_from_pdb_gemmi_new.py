import gemmi
from pathlib import Path
from tqdm import tqdm
from foldeverything.data import const

def extract_seqs(pdb_file: str) -> None:
    """Extract sequences from a PDB file."""
    # pdb_file: Path = Path(pdb_file)
    # pdb_id = pdb_file.stem.split(".")[0].lower()

    structure = gemmi.read_structure(pdb_file)
    # structure.setup_entities()
    # structure.add_entity_types()
    # structure.assign_subchains()
    # structure.ensure_entities()

    structure.merge_chain_parts()
    structure.remove_waters()
    structure.remove_hydrogens()
    structure.remove_alternative_conformations()
    structure.remove_empty_chains()
    try:
        structure.transform_to_assembly(
            structure.assemblies[0].name,
            gemmi.HowToNameCopiedChain.AddNumber,
        )
    except Exception as e:
        print(f"Warning: fail to transform_to_assembly: {pdb_file}")

    # prots, rnas, dnas, ligands = [], [], [], []
    prots = dict()
    prots_pdbx = dict()
    for entity_id, entity in enumerate(structure.entities):
        # print(entity_id, entity, gemmi.EntityType.Polymer, entity.polymer_type.name)
        # exit()
        if "1zea" in pdb_file:
            print(entity_id, entity, gemmi.EntityType.Polymer, entity.polymer_type.name)
        entity: gemmi.Entity
        # print(entity_id)
        # print(gemmi.EntityType.Polymer, entity.entity_type, entity.polymer_type.name)

        # if entity.entity_type == gemmi.EntityType.Polymer and (
        #     entity.polymer_type.name in {"PeptideL", "Rna", "Dna"}
        # ):
        if entity.entity_type == gemmi.EntityType.Polymer:
            # Fetch the sequence
            seq = entity.full_sequence
            seq = [gemmi.Entity.first_mon(item) for item in seq]
            
            # if entity.polymer_type.name == "PeptideL":
            
            seq_pdbx = gemmi.pdbx_one_letter_code(seq, gemmi.ResidueKind.AA)
            prots_pdbx[entity.name] = seq_pdbx
            
            seq_boltz = []
            for aa in seq:
                aa_name_boltz_extended = const.protein_letters_3to1_extended.get(aa, "X")
                if aa_name_boltz_extended == "X":
                    aa_name = gemmi.find_tabulated_residue(aa).one_letter_code
                else:
                    aa_name = aa_name_boltz_extended
                if not aa_name.strip():
                    seq_boltz.append("X")
                else:
                    seq_boltz.append(aa_name)
            seq_boltz = "".join(seq_boltz)
            prots[entity.name] = seq_boltz
            
    return prots, prots_pdbx

if __name__ == "__main__":
    
    pdb_dir = Path("/data/scratch/wxsh/data/antibodies/sabdab/all_structures/raw")
    output_fn = "sequences_from_raw_pdb_boltz_extended.fasta"
    output_fn_pdbx = "sequences_from_raw_pdb_pdbx_code.fasta"

    pdb_files = list(pdb_dir.rglob("*.pdb"))
    with open(output_fn, "w") as fout, open(output_fn_pdbx, "w") as fout_pdbx:
        for file in tqdm(pdb_files):
            pdb_id = file.stem.lower()
            sequences, sequences_pdbx = extract_seqs(str(file))
            for chain in sequences:
                if sequences[chain]:
                    fout.write(f">{pdb_id.lower()}_{chain}\n{sequences[chain]}\n")
            for chain in sequences_pdbx:
                if sequences_pdbx[chain]:
                    fout_pdbx.write(f">{pdb_id.lower()}_{chain}\n{sequences_pdbx[chain]}\n")
