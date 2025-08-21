"""Create a mapping from structure and chain ID to MSA indices."""

import argparse
import json
import multiprocessing
from pathlib import Path

import gemmi
from tqdm import tqdm

from foldeverything.data import const


def get_release_date(block: gemmi.cif.Block) -> str:
    """Get the released date.

    Parameters
    ----------
    block : gemmi.cif.Block
        The block to process.

    Returns
    -------
    str
        The released date.

    """
    revision = "_pdbx_audit_revision_history.revision_date"
    release_date = block.find([revision])[0][0]
    return release_date


def extract_seqs(pdb_file: str) -> None:
    """Extract sequences from a PDB file."""
    pdb_file: Path = Path(pdb_file)
    pdb_id = pdb_file.stem.split(".")[0].lower()
    try:
        # Parse MMCIF input file
        block = gemmi.cif.read(str(pdb_file))[0]
        release_date = get_release_date(block)
        structure = gemmi.make_structure_from_block(block)

    except RuntimeError:
        print(f"Failed to read {pdb_file}")
        return {"id": pdb_id, "proteins": [], "rnas": [], "dnas": [], "ligands": []}

    structure.merge_chain_parts()
    structure.remove_waters()
    structure.remove_hydrogens()
    structure.remove_alternative_conformations()
    structure.remove_empty_chains()
    structure.transform_to_assembly(
        structure.assemblies[0].name,
        gemmi.HowToNameCopiedChain.AddNumber,
    )

    prots, rnas, dnas, ligands = [], [], [], []
    for entity_id, entity in enumerate(structure.entities):
        entity: gemmi.Entity
        if entity.entity_type == gemmi.EntityType.Polymer and (
            entity.polymer_type.name in {"PeptideL", "Rna", "Dna"}
        ):
            # Fetch the sequence
            seq = entity.full_sequence
            seq = [gemmi.Entity.first_mon(item) for item in seq]
            subchains = ",".join(entity.subchains)

            if entity.polymer_type.name == "PeptideL":
                seq_boltz_old = [const.prot_token_to_letter.get(aa, "X") for aa in seq]
                seq_boltz_old = "".join(seq_boltz_old)

                seq_boltz_new = [const.protein_letters_3to1_extended.get(aa, "X") for aa in seq]
                seq_boltz_new = "".join(seq_boltz_new)

                residue_aas = []
                residue_aas_boltz_extended_and_gemmi = []
                for aa in seq:
                    ### merge knowledge from boltz and Gemmi:
                    aa_name_boltz_extended = const.protein_letters_3to1_extended.get(aa, "X")
                    if aa_name_boltz_extended == "X":
                        aa_name = gemmi.find_tabulated_residue(aa).one_letter_code
                    else:
                        aa_name = aa_name_boltz_extended
                    if not aa_name.strip():
                        residue_aas_boltz_extended_and_gemmi.append("X")
                    else:
                        residue_aas_boltz_extended_and_gemmi.append(aa_name)

                    ### ori:
                    aa_name = gemmi.find_tabulated_residue(aa).one_letter_code
                    if not aa_name.strip():
                        residue_aas.append("X")
                    else:
                        residue_aas.append(aa_name)
                seq_gemmi = "".join(residue_aas)
                seq_boltz_extended_and_gemmi = "".join(residue_aas_boltz_extended_and_gemmi)
                seq_pdbx = gemmi.pdbx_one_letter_code(seq, gemmi.ResidueKind.AA)

                # print("seq_boltz_old", seq_boltz_old)
                # print("seq_boltz_new", seq_boltz_new)
                # print("seq_gemmi", seq_gemmi)
                # print("seq_pdbx", seq_pdbx)
                # print(seq_boltz_extended_and_gemmi)
                # exit()

                # prots.append((pdb_id, entity_id, entity.name, subchains, seq))
                prots.append((pdb_id, entity_id, entity.name, subchains, \
                    {
                        "seq_boltz": seq_boltz_old,
                        "seq_boltz_extended": seq_boltz_new,
                        "seq_gemmi": seq_gemmi,
                        "seq_pdbx": seq_pdbx,
                        "seq_boltz_extended_and_gemmi": seq_boltz_extended_and_gemmi,
                    }))
            elif entity.polymer_type.name == "Rna":
                seq = [const.rna_token_to_letter.get(aa, "N") for aa in seq]
                seq = "".join(seq)
                rnas.append((pdb_id, entity_id, entity.name, subchains, seq))
            elif entity.polymer_type.name == "Dna":
                seq = [const.dna_token_to_letter.get(aa, "N") for aa in seq]
                seq = "".join(seq)
                dnas.append((pdb_id, entity_id, entity.name, subchains, seq))

        elif entity.entity_type in (
            gemmi.EntityType.NonPolymer,
            gemmi.EntityType.Branched,
        ):
            # Fetch the residues from the first subchain, if any
            if not entity.subchains:
                continue
            subchain_id = entity.subchains[0]
            subchain = next(
                c for c in structure[0].subchains() if c.subchain_id() == subchain_id
            )
            # Create a fake sequence by concatenating the residue names
            seq = "_".join([res.name for res in subchain])
            subchains = ",".join(entity.subchains)
            ligands.append((pdb_id, entity_id, entity.name, subchains, seq))

    return {
        "id": pdb_id,
        "released": release_date,
        "proteins": prots,
        "rnas": rnas,
        "dnas": dnas,
        "ligands": ligands,
    }


def main(
    pdb_dir: Path,
    outdir: Path,
    num_processes: int,
    subfolders: bool = True,
) -> None:
    """Create mapping."""
    # Get all PDB files
    print("Looking for PDB files...")
    if subfolders:
        folders = [item for item in pdb_dir.iterdir() if item.is_dir()]
    else:
        folders = [pdb_dir]
    files = [str(item) for f in folders for item in f.glob("*.cif.gz")]

    print("len(files):", len(files))

    # Extract sequences
    all_data = []
    num_processes = min(num_processes, multiprocessing.cpu_count())

    if num_processes == 1:
        for pdb_file in tqdm(files):
            data = extract_seqs(pdb_file)
            all_data.append(data)
    else:
        with multiprocessing.Pool(num_processes) as pool:  # noqa: SIM117
            with tqdm(total=len(files)) as pbar:
                for data in pool.imap_unordered(extract_seqs, files):
                    all_data.append(data)
                    pbar.update()

    # Write sequences to fasta files
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "sequences.json").open("w") as f:
        json.dump(all_data, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--num_processes", type=int, default=multiprocessing.cpu_count()
    )
    parser.add_argument("--subfolders", action="store_true")
    args = parser.parse_args()
    main(**vars(args))