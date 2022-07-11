import pytest
from src.classify_homologues.utils import *


def test_parser() -> None:
    """
    Test to check arguments are correctly specified.
    """
    assert 1 == 1

def test_read_smiles_csv() -> None:
    """
    Test to check there are enough molecules to work with for homologous series classification (minimum 2).
    :return: None
    """
    smiles, mols, smiles_torem,idxtorem = read_smiles_csv("tests/test1_smiles_23.csv")
    len_mols = len(mols)
    assert len_mols >= 2 , "Not enough molecules. Check validity of SMILES input or add more SMILES."


def test_read_labels_csv() -> None:
    """
    Test to check every molecule has a name i.e. one label per Mol.
    """
    idxtorem = []
    labels = read_labels_csv("tests/test1_labels_23.csv", idxtorem)
    assert len(labels) == 23, "Too many/few labels. Compare lengths of SMILES and labels inputs."





#def test_setup_repeating_unit() -> None:
#    """
#    Test to check that repeating units were built correctly, end-to-end.
#    """
