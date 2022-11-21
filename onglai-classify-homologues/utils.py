from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.rdmolops
import numpy as np
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
from itertools import compress, repeat
import datamol as dm
import pandas as pd

try:
    import Image
except ImportError:
    from PIL import Image


def read_input_csv_smiles_name(path_to_csv: str, smiles_col="SMILES", name_col="Name"):
    """Function to parse input CSV, containing minimum 2 columns containing SMILES and Name."""

    # check column inputs - str or int
    # if smiles_col

    cols = [smiles_col, name_col]
    input_df = pd.read_csv(path_to_csv, usecols=cols)

    # parse SMILES
    smiles = input_df[smiles_col]
    mols = [Chem.MolFromSmiles(s) for s in smiles]

    # check validity of SMILES - save invalid SMILES idxs
    idxtorem = [i for i, j in enumerate(mols) if j is None]  # indexes of empty mols
    smiles_torem = list()
    if len(idxtorem) > 0:
        print(
            str(len(idxtorem))
            + " parsed SMILES gave empty Mol objects (check validity of SMILES!). Removing those Mols."
        )
        smiles_torem = [y for x, y in enumerate(smiles) if x in idxtorem]
        smiles = [y for x, y in enumerate(smiles) if x not in idxtorem]
        mols = [y for x, y in enumerate(mols) if x not in idxtorem]

    # parse names
    labels = input_df[name_col].to_list()
    labels = [y for x, y in enumerate(labels) if x not in idxtorem]

    # generate df for generate_df
    df = pd.DataFrame({"SMILES": smiles, "Mols": mols, "Labels": labels})

    return smiles, mols, smiles_torem, idxtorem, labels, input_df, df, path_to_csv


# def read_smiles_csv(path_to_smiles_csv: str): #sys.argv[1]
#     '''Function to parse SMILES from CSV file and create molecule objects. Isolates unparseable SMILES and returns them, and their idxs as a list.'''
#     with open(path_to_smiles_csv) as f:
#         smiles = [line.strip().replace('"','') for line in f]
#         mols = [AllChem.MolFromSmiles(smile) for smile in smiles]
#         #check validity of SMILES - if there are any empty mol objects caused by unparseable SMILES
#         idxtorem = [i for i,j in enumerate(mols) if j is None] #indexes of empty mols
#         smiles_torem = list()
#         if len(idxtorem) > 0:
#             print(str(len(idxtorem))+ " parsed SMILES gave empty Mol objects (check validity of SMILES!). Removing those Mols.")
#             smiles_torem = [y for x,y in enumerate(smiles) if x in idxtorem]
#             smiles = [y for x,y in enumerate(smiles) if x not in idxtorem]
#             mols = [y for x,y in enumerate(mols) if x not in idxtorem]
#         return smiles, mols, smiles_torem, idxtorem


def write_removed_smiles(smiles_torem):
    with open("output/removed_smiles.txt", "w") as text_file:
        text_file.write("\n".join(smiles_torem))


# def read_labels_csv(path_to_labels_csv, idxtorem): #sys.argv[2]
#     '''Function to read in list of labels corresponding to SMILES.'''
#     with open(path_to_labels_csv) as f:
#         labels = [line.strip().replace('"','') for line in f]
#         labels = [j for i,j in enumerate(labels) if i not in idxtorem]
#     return labels


def setup_repeating_unit(smarts, min, max):
    """Function to generate list of RU chains as query molecules from SMARTS."""
    smiles_ru = []
    [smiles_ru.append(smarts * i) for i in range(min, max + 1)]  # add min
    smiles_ru = [x[:-1] for x in smiles_ru]  # remove last hyphen in each string
    ru = [Chem.MolFromSmarts(smi) for smi in smiles_ru]
    return ru
    # return(smiles_ru)


def detect_repeating_units(mols, labels, ru):
    """Function to detect whether molecules contain repeating units."""
    # set up RU-match matrix for detection of RU in mols
    mat1 = SubstructMatchMatrix_ru_mols(mols, ru, accountForRings=True)
    mat_array_sums = []
    for x, y in enumerate(mat1):
        mat_array_sums.append(int(np.sum(mat1[x])))
    n_mols_no_ru = mat_array_sums.count(0)
    if n_mols_no_ru > 0:
        print(
            str(n_mols_no_ru)
            + " mols have no repeating unit chains of minimum length specified."
        )
    # remove mols with no RU matches
    fil_ru = [
        bool(x) for x in mat_array_sums
    ]  # those which are False have array_sum = 0 i.e. no alkyls
    mols_no_ru_matches = list(compress(mols, [not i for i in fil_ru]))
    labels_mols_no_ru_matches = list(compress(labels, [not i for i in fil_ru]))
    # if len(mols_no_ru_matches) > 0:
    #    nans = DrawMolsZoomed(mols=mols_no_ru_matches, legends=labels_mols_no_ru_matches, molsPerRow=5)
    #    nans.save("output_rmdum_tmf/no_min_ru_matches.png")
    mols_with_ru = list(compress(mols, fil_ru))
    labels_mols_with_ru = list(compress(labels, fil_ru))
    return (
        mols_no_ru_matches,
        labels_mols_no_ru_matches,
        mols_with_ru,
        labels_mols_with_ru,
    )


def hasSubstructMatchAccountForRings(mol, q):
    """Function to exclude substructure matches if they are part of rings. Credit Paolo Tosco."""
    matches = mol.GetSubstructMatches(q)
    hasMatch = False
    for match in matches:
        hasMatch = True
        for i, j in enumerate(match):
            if q.GetAtomWithIdx(i).IsInRing() ^ mol.GetAtomWithIdx(j).IsInRing():
                hasMatch = False
                break
        if hasMatch:
            break
    return hasMatch


# from https://gist.github.com/ptosco/26af473fc1f3129878ca86cb070afe3a


def SubstructMatchMatrix_ru_mols(mols, ru, accountForRings=True):
    """Function to generate matrix of 1s and 0s for RU substructure matching in molecules."""
    mat = np.zeros((len(mols), len(ru)))
    for a, b in enumerate(mols):
        for i, j in enumerate(ru):
            if accountForRings:  # exclude RU matches if they are part of rings
                mat[a, i] = hasSubstructMatchAccountForRings(b, j)
            else:
                mat[a, i] = mols[a].HasSubstructMatch(j)
    return mat


def replacecore_detect_homologue_cores(mols_with_ru, ru):
    """Function that performs fragmentation (RU matching and replacement with a dummy atom) to isolate core fragments in molecule object. Dummies are removed."""
    mat2 = SubstructMatchMatrix_ru_mols(
        mols_with_ru, ru, accountForRings=True
    )  # set up RU-match matrix
    patts, cores = replacecore_longest_RU_match(mols_with_ru, mat2, ru)  #
    patts = [dm.remove_dummies(m, dummy="*") for m in patts]  # remove dummies
    cores = [dm.remove_dummies(n, dummy="*") for n in cores]  # remove dummies
    return patts, cores


def replacecore_keepdummies_detect_homologue_cores(mols_with_ru, ru):
    """Function that performs fragmentation (RU matching and replacement with a dummy atom) to isolate core fragments in molecule object. Dummies NOT removed."""
    mat2 = SubstructMatchMatrix_ru_mols(
        mols_with_ru, ru, accountForRings=True
    )  # set up RU-match matrix
    patts, cores = replacecore_longest_RU_match(mols_with_ru, mat2, ru)  #
    return patts, cores


def fragment_into_cores(mols_with_ru, ru, frag_steps):
    """Function that repeats replacecore_detect_homologue_cores n times, taking output of last as input of next fragmentation step. List of empty_cores_idx generated at the end of all steps."""
    lists_patts = []
    lists_cores = []
    for i in range(frag_steps):
        if i == 0:
            m = mols_with_ru
        else:
            m = lists_cores[-1]  # take output of previous iteration as input
        patts, cores = replacecore_detect_homologue_cores(m, ru)
        lists_patts.append(patts)
        lists_cores.append(cores)
    # if there are any empty cores after all frag steps, get their index. Else, empty_cores_idx is itself empty list.
    empty_cores_idx = [i for i, j in enumerate(lists_cores[-1]) if j.GetNumAtoms() == 0]
    return lists_patts, lists_cores, empty_cores_idx


def my_ReplaceCore(mol, sidechain):
    """Function that returns the original input mol if there is no matching sidechain present."""
    rc = Chem.ReplaceCore(mol, sidechain)
    if rc is None:  # mol does not contain sidechain whatsoever
        return mol
    else:
        return rc


def replacecore_longest_RU_match(mols, mat, ru):
    """Function to replace the longest RU substructure match from each molecule with dummy atoms using RDKit's Chemical Transformations functionality. Returns remaining cores and the RU replaced by a dummy atom."""
    cores = list()
    patts = list()
    for x, y in enumerate(mols):
        patts.append(ru[int(np.sum(mat[x]) - 1)])  # zero-index
        cores.append(my_ReplaceCore(y, patts[x]))
    return patts, cores


def detect_mols_made_of_ru(mols_with_ru, labels_mols_with_ru, empty_cores_idx):
    """Function to detect, depict, then discard molecules made solely of RUs such as PEGs. Depictions written as png outputs."""
    mols_made_of_ru = [
        j for i, j in enumerate(mols_with_ru) if (i in empty_cores_idx)
    ]  # isolate Mol and Label with empty core after first chopping i.e. entire mol is made of ru
    labels_made_of_ru = [
        j for i, j in enumerate(labels_mols_with_ru) if (i in empty_cores_idx)
    ]
    mols_to_classify = [
        j for i, j in enumerate(mols_with_ru) if (i not in empty_cores_idx)
    ]
    labels_to_classify = [
        j for i, j in enumerate(labels_mols_with_ru) if (i not in empty_cores_idx)
    ]
    if len(mols_made_of_ru) > 0:
        #    pure_repeating_units = DrawMolsZoomed(mols_made_of_ru, labels_made_of_ru)
        #    pure_repeating_units.save("output_rmdum_tmf/mols_pure_repeating_units.png")
        print(
            str(len(mols_made_of_ru))
            + " molecule(s) are made purely of repeating units of minimum length specified (default=3)."
        )
    else:
        print("No mols made of pure RUs.")
    return (
        mols_made_of_ru,
        labels_made_of_ru,
        mols_to_classify,
        labels_to_classify,
    )  # includes non-series containing RU


def process_patts_cores(
    lists_patts, lists_cores, empty_cores_idx
):  # prepare dataframe???
    """Filter out empty cores from patterns and cores lists in preparation for df assembly."""
    for i, j in enumerate(lists_patts):  # lists_patts should be same len as lists_cores
        lists_patts[i] = [
            q for p, q in enumerate(lists_patts[i]) if (p not in empty_cores_idx)
        ]
        # for i,j in enumerate(lists_cores):
        lists_cores[i] = [
            q for p, q in enumerate(lists_cores[i]) if (p not in empty_cores_idx)
        ]
    return lists_patts, lists_cores


def generate_df(
    lists_patts, lists_cores, mols_to_classify, labels_to_classify, df, mols_made_of_ru
):
    canosmiles_finalcores = [
        Chem.MolToSmiles(c, canonical=True) for c in lists_cores[-1]
    ]
    nonempty_df = pd.DataFrame(
        {
            "Mols": mols_to_classify,
            "Labels": labels_to_classify,
            "CanoSmiles_FinalCores": canosmiles_finalcores,
        }
    )
    # assign series numbers 0 to n
    nonempty_df["SeriesNo"] = (
        nonempty_df.groupby("CanoSmiles_FinalCores")
        .filter(lambda group: len(group) > 1)
        .groupby("CanoSmiles_FinalCores")
        .ngroup()
    )
    result_df = pd.merge(
        df, nonempty_df, how="left", on=["Mols", "Labels"]
    )  # add back mols previously filtered out: mols_no_ru_matches & mols_made_of_ru
    result_df["SeriesNo"] = result_df["SeriesNo"].fillna(
        -1
    )  # encode -1 for mols_no_ru_matches & mols_made_of_ru & nonseries with RU matches
    result_df = encode_pureRU_minus2(
        result_df, mols_made_of_ru
    )  # encode -2 for mols_made_of_ru
    result_df["SeriesNo"] = np.where(
        (result_df["SeriesNo"] == -1) & (result_df["CanoSmiles_FinalCores"].notnull()),
        -3,
        result_df["SeriesNo"],
    )  # encode -3 for non-series with RU
    result_df.SeriesNo = result_df.SeriesNo.astype(int)
    classified_series = result_df[result_df["SeriesNo"] > -1]
    return classified_series, result_df


def encode_pureRU_minus2(result_df, mols_made_of_ru):
    rdkitsmiles_mols_made_of_ru = [
        Chem.MolToSmiles(m, canonical=True) for m in mols_made_of_ru
    ]
    rdkitsmiles_mols = [
        Chem.MolToSmiles(m, canonical=True) for m in result_df.Mols
    ]  # these are different from input SMILES!
    idx_pureRU = [
        i for i, j in enumerate(rdkitsmiles_mols) if j in rdkitsmiles_mols_made_of_ru
    ]
    result_df["SeriesNo"] = np.where(
        result_df.index.isin(idx_pureRU), -2, result_df["SeriesNo"]
    )
    return result_df


def detect_mols_nonseries(result_df):
    # encode one_member_series with -2 in SeriesNo
    nonseries = result_df.loc[
        (result_df["SeriesNo"] == -3) & (result_df["CanoSmiles_FinalCores"].notnull())
    ]
    if len(nonseries.Mols) > 0:
        mols_nonseries = [i for i in nonseries.Mols]
        labs_nonseries = [i for i in nonseries.Labels]
        # pl_onememseries = DrawMolsZoomed(mols_onememseries,labs_onememseries,molsPerRow=5)
        # pl_onememseries.save("output_rmdum_tmf/onememseries_containing_repeating_unit.png")
        print(
            str(len(nonseries.Mols))
            + " molecule(s) contain RUs of minimum chain length specified but do not form series."
        )
    else:
        mols_nonseries = []
        labs_nonseries = []
        print(
            str(len(nonseries.Mols))
            + " molecule(s) contain RUs of minimum chain length specified but have unique cores (one-member series)."
        )

    return mols_nonseries, labs_nonseries, nonseries


def detect_cores_classified_series(classified_series):
    grpdmols = classified_series.groupby("CanoSmiles_FinalCores").SMILES.apply(
        list
    )  # lists of SMILES of molecules in each series
    for i, j in enumerate(grpdmols):
        grpdmols[i] = grpdmols[i] + [grpdmols.keys()[i]]
    return grpdmols


def depict_cores_summary(grpdmols):
    """Depict cores of classified series."""
    if len(grpdmols.keys()) > 0:
        final_cores = [Chem.MolFromSmiles(i) for i in grpdmols.keys()]
        leg_final_cores = [str(idx) for idx, y in enumerate(grpdmols.keys())]
        # cores_summary = DrawMolsZoomed(final_cores, legends=leg_final_cores, molsPerRow=5)
        # cores_summary.save("output_rmdum_tmf/cores_summary.png")
    else:
        final_cores = []
        leg_final_cores = []

    return final_cores, leg_final_cores


def generate_classified_series_summary(result_df):
    """Write classified series results as CSV."""
    inchis = [Chem.inchi.MolToInchi(i) for i in result_df.Mols]
    inchikeys = [Chem.inchi.MolToInchiKey(i) for i in result_df.Mols]
    mf = [Chem.rdMolDescriptors.CalcMolFormula(i) for i in result_df.Mols]
    monoiso_mass = [round(Chem.Descriptors.ExactMolWt(i), 4) for i in result_df.Mols]
    out = result_df[["SeriesNo", "Labels", "CanoSmiles_FinalCores", "SMILES",]].copy()
    out["InChI"] = inchis
    out["InChIKey"] = inchikeys
    out["molecular_formula"] = mf
    out["monoisotopic_mass"] = monoiso_mass
    out.rename(
        columns={
            "SeriesNo": "series_no",
            "Labels": "cpd_name",
            "canoSMILES_molfrags": "core_fragments",
        },
        inplace=True,
    )
    out.to_csv("output/" + "classified_series.csv", index=False)


def depict_classified_series(grpdmols, classified_series):
    grpdmols = [[Chem.MolFromSmiles(s) for s in g] for g in grpdmols]
    lgs = [
        i for i in classified_series.groupby("CanoSmiles_FinalCores").Labels.apply(list)
    ]
    for i, j in enumerate(lgs):
        lgs[i] = lgs[i] + ["core"]
    list_grid_images = []
    for i, j in enumerate(grpdmols):
        list_grid_images.append(
            DrawMolsZoomed(grpdmols[i], legends=lgs[i], molsPerRow=5)
        )
    # save each plot per group
    [img.save("output" + str(idx) + ".png") for idx, img in enumerate(list_grid_images)]


def DrawMolsZoomed(mols, legends, molsPerRow=3, subImgSize=(300, 300)):  # , leg):
    """Function to draw rows of zoomed molecules. Credit Rocco Moretti."""
    nRows = len(mols) // molsPerRow
    if len(mols) % molsPerRow:
        nRows += 1
    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
    full_image = Image.new("RGB", fullSize)
    for ii, mol in enumerate(mols):
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)
        le = legends[ii]
        column = ii % molsPerRow
        row = ii // molsPerRow
        offset = (column * subImgSize[0], row * subImgSize[1])
        d2d = rdMolDraw2D.MolDraw2DCairo(subImgSize[0], subImgSize[1])
        d2d.DrawMolecule(mol, legend=le)
        d2d.FinishDrawing()
        sub = Image.open(BytesIO(d2d.GetDrawingText()))
        full_image.paste(sub, box=offset)
    return full_image


# from https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/CAHGTkV8sdfb4Q7FLn9C5MTwrqiJjHtnXK%2Bmz2SY3_4j2eAtevQ%40mail.gmail.com/#msg36477772


def delete_longest_RU_match(mols, mat, ru):
    """Function to delete the longest RU substructure match from each molecule. Returns remaining cores and the RU deleted."""
    cores = list()
    patt = list()
    for x, y in enumerate(mols):
        patt.append(ru[int(np.sum(mat[x]) - 1)])
        cores.append(AllChem.DeleteSubstructs(y, patt[x]))
    return patt, cores


def replace_longest_RU_match(mols, mat, ru):
    """Function to replace the longest RU substructure match from each molecule with C. Returns remaining cores and the RU replaced by C."""
    cores = list()
    patt = list()
    for x, y in enumerate(mols):
        patt.append(ru[int(np.sum(mat[x]) - 1)])
        cores.append(
            AllChem.ReplaceSubstructs(
                y, patt[x], Chem.MolFromSmiles("C"), replaceAll=True
            )[0]
        )
    return patt, cores


def GetNumFrags(mol):
    return len(Chem.GetMolFrags(mol, asMols=True))


def print_output_summary(result_df, nonseries, mols_no_ru_matches, mols_made_of_ru):
    num_series = result_df.SeriesNo.max() + 1  # because zero-indexed
    if num_series < 0:
        num_series = 0
    mols_classified = (
        len(result_df.Mols)
        - len(nonseries.Mols)
        - len(mols_no_ru_matches)
        - len(mols_made_of_ru)
    )
    return num_series, mols_classified


def generate_output_summary(
    mols_classified,
    num_series,
    ru_in,
    mols_no_ru_matches,
    nonseries,
    mols_made_of_ru,
    min_length,
    max_length,
    frag_steps,
    runtime,
    path_to_csv,
):
    ru_in = ru_in[:-1]
    mols_no_ru_matches = len(mols_no_ru_matches)
    nonseries = len(nonseries)
    mols_made_of_ru = len(mols_made_of_ru)
    with open("output/classification-results.txt", "w") as text_file:
        text_file.write(
            "Homologue classification for %s complete! \nRepeating Unit (RU) = %s.\nFragmentation steps = %s.\n%s molecules were classified into %s homologous series. [SeriesNo >= 0]\n"
            % (path_to_csv, ru_in, frag_steps, mols_classified, num_series)
        )
        text_file.write(
            "Molecules with no RU matches of minimum length %s units, maximum length %s units: %s. [SeriesNo = -1]\n"
            % (min_length, max_length, mols_no_ru_matches)
        )
        text_file.write(
            "Molecules consisting of purely RUs: %s. [SeriesNo = -2]\n"
            % (mols_made_of_ru)
        )
        text_file.write(
            "Molecules having RU matches of minimum length %s units but do not form series: %s. [SeriesNo = -3]\n"
            % (min_length, nonseries)
        )
        text_file.write("Algorithm runtime in seconds: %s.\n" % (runtime))
