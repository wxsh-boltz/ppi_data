from dataclasses import dataclass
from typing import List, Dict

@dataclass
class Dataset:
    dataset_name: str
    # dataset_type: str
    json_path: str
    alignment_paths: Dict[str, List[str]] | None = None
    sequence_identity_threshold: Dict[str, float] | float = None
    edit_distance_threshold: Dict[str, float] | float = None

@dataclass
class GeneralPPIDataset(Dataset):
    # alignment_paths: List[str]
    seq_path: str = None
    def __post_init__(self):
        self.dataset_type = "general_ppi"

@dataclass
class AntibodyDataset(Dataset):
    # target_to_rcsb_train_mmseqs_alignment_paths: List[str] | None = None
    # target_to_target_mmseqs_alignment_path: str | None = None
    target_peptide_seq_path: str | None = None
    # target_peptide_to_rcsb_train_edit_distance_path: str | None = None
    target_seq_path: str = None
    vl_seq_path: str = None
    vh_seq_path: str = None

    def __post_init__(self):
        self.dataset_type = "antibody"
        self.target_role = "antigen"        


@dataclass
class MiniproteinDataset(Dataset):
    target_seq_path: str = None
    binder_seq_path: str = None
    target_peptide_seq_path: str | None = None

    def __post_init__(self):
        self.dataset_type = "miniprotein"
        self.target_role = "target"   


@dataclass
class TCRpMHCDataset(Dataset):
    mhc_seq_path: str = None
    tcra_seq_path: str = None
    tcrb_seq_path: str = None
    peptide_seq_path: str = None

    def __post_init__(self):
        self.dataset_type = "tcr_phmc"
        self.target_role = "peptide"  