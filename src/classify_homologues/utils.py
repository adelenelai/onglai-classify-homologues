from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.rdmolops
import numpy as np
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
from itertools import compress, repeat
import datamol as dm
try:
    import Image
except ImportError:
    from PIL import Image

def read_smiles_csv(path_to_smiles_csv): #sys.argv[1]
    '''Function to read in list of SMILES and create molecule objects.'''
    with open(path_to_smiles_csv) as f:
        smiles = [line.strip().replace('"','') for line in f]
        mols = [AllChem.MolFromSmiles(smile) for smile in smiles]
    return smiles, mols

def read_labels_csv(path_to_labels_csv): #sys.argv[2]
    '''Function to read in list of labels corresponding to SMILES.'''
    with open(path_to_labels_csv) as f:
        labels = [line.strip().replace('"','') for line in f]
    return labels

def setup_repeating_unit(smarts, min, max):
    '''Function to generate list of RU chains as query molecules from SMARTS.'''
    smiles_ru = []
    [smiles_ru.append(smarts*i) for i in range(min, max+1)] #add min
    smiles_ru = [x[:-1] for x in smiles_ru] #remove last hyphen in each string
    ru = [Chem.MolFromSmarts(smi) for smi in smiles_ru]
    return(ru)
    #return(smiles_ru)


def detect_repeating_units(mols, labels, ru):
    '''Function to detect whether molecules contain repeating units.'''
    # set up RU-match matrix for detection of RU in mols
    mat1 = SubstructMatchMatrix_ru_mols(mols, ru, accountForRings=True)
    mat_array_sums = []
    length_ru_chains_to_chop = []
    for x, y in enumerate(mat1):
        mat_array_sums.append(int(np.sum(mat1[x])))
        length_ru_chains_to_chop.append(mat_array_sums[x] + 2)  # n = sum + 2, where n is the length of C chain to chop
    n_mols_no_ru = mat_array_sums.count(0)
    if n_mols_no_ru>0:
        print(str(n_mols_no_ru) + " mols have repeating unit chains but are below minimum length specified (default = 3).")
    #remove mols with no RU matches
    fil_ru = [bool(x) for x in mat_array_sums] #those which are False have array_sum = 0 i.e. no alkyls
    mols_no_ru_matches = list(compress(mols, [not i for i in fil_ru]))
    labels_mols_no_ru_matches = list(compress(labels, [not i for i in fil_ru]))
    if len(mols_no_ru_matches) > 0:
        nans = DrawMolsZoomed(mols=mols_no_ru_matches, legends=labels_mols_no_ru_matches, molsPerRow=5)
        nans.save("output_rmdum_tmf/no_repeating_unit_matches.png")
    mols_with_ru = list(compress(mols, fil_ru))
    labels_mols_with_ru = list(compress(labels, fil_ru))
    return mols_no_ru_matches, labels_mols_no_ru_matches, mols_with_ru, labels_mols_with_ru

# def detect_homologue_cores(mols_with_ru, ru):
#     '''Function that performs RU matching-and-removal twice to isolate/detect cores in molecule object. Idxs of empty cores generated.'''
#     mat2 = SubstructMatchMatrix_ru_mols(mols_with_ru, ru, accountForRings=True) #set up RU-match matrix for 1st RU removal from mols with RU
#     patt1, cores1 = delete_longest_RU_match(mols_with_ru, mat2, ru) #first removal
#     empty_cores_idx = [i for i, j in enumerate(cores1) if j.GetNumAtoms() == 0] #isolate empty cores' idxs after first chopping, occur when mol is 100% made of RU
#     mat3 = SubstructMatchMatrix_ru_mols(cores1, ru, accountForRings=True) #set up RU-match matrix for 2nd RU removal from cores1
#     patt2, cores2 = delete_longest_RU_match(cores1, mat3, ru) #second removal
#     return patt1, cores1, patt2, cores2, empty_cores_idx


# def replaceRU_detect_homologue_cores(mols_with_ru, ru):
#     '''Function that performs fragmentation (RU matching and replacement with a dummy atom) n times (n = frag_steps) to isolate core fragments in molecule object. Idxs of empty cores generated at the end.'''
#     mat2 = SubstructMatchMatrix_ru_mols(mols_with_ru, ru, accountForRings=True) #set up RU-match matrix for 1st RU removal from mols with RU
#     patt1, cores1 = replace_longest_RU_match(mols_with_ru, mat2, ru) #first removal
#     empty_cores_idx = [i for i, j in enumerate(cores1) if j.GetNumAtoms() == 0] #isolate empty cores' idxs after first chopping, occur when mol is 100% made of RU
#     mat3 = SubstructMatchMatrix_ru_mols(cores1, ru, accountForRings=True) #set up RU-match matrix for 2nd RU removal from cores1
#     patt2, cores2 = replace_longest_RU_match(cores1, mat3, ru) #second removal
#     return patt1, cores1, patt2, cores2, empty_cores_idx

def hasSubstructMatchAccountForRings(mol, q):
    """Function to exclude substructure matches if they are part of rings. Credit Paolo Tosco."""
    matches = mol.GetSubstructMatches(q)
    hasMatch = False
    for match in matches:
        hasMatch = True
        for i, j in enumerate(match):
            if (q.GetAtomWithIdx(i).IsInRing() ^ mol.GetAtomWithIdx(j).IsInRing()):
                hasMatch = False
                break
        if (hasMatch):
            break
    return hasMatch
#from https://gist.github.com/ptosco/26af473fc1f3129878ca86cb070afe3a


def SubstructMatchMatrix_ru_mols(mols, ru, accountForRings=True):
    '''Function to generate matrix of 1s and 0s for RU substructure matching in molecules.'''
    mat = np.zeros((len(mols), len(ru)))
    for a,b in enumerate(mols):
        for i,j in enumerate(ru):
            if accountForRings:
                mat[a,i] = hasSubstructMatchAccountForRings(b,j)
            else:
                mat[a,i] = mols[a].HasSubstructMatch(j)
    return mat

def replacecore_detect_homologue_cores(mols_with_ru, ru):
    '''Function that performs fragmentation (RU matching and replacement with a dummy atom) to isolate core fragments in molecule object.'''
    mat2 = SubstructMatchMatrix_ru_mols(mols_with_ru, ru, accountForRings=True) #set up RU-match matrix
    patts, cores = replacecore_longest_RU_match(mols_with_ru, mat2, ru) #
    patts = [dm.remove_dummies(m,dummy='*') for m in patts] #remove dummies
    cores = [dm.remove_dummies(n,dummy='*') for n in cores] #remove dummies
    return patts, cores

def fragment_into_cores(mols_with_ru, ru, frag_steps):
    '''Function that repeats replacecore_detect_homologue_cores n times, taking output of last as input of next fragmentation step. List of empty_cores_idx generated at the end of all steps.'''
    lists_patts = []
    lists_cores = []
    for i in range(frag_steps):
        if i == 0:
            m = mols_with_ru
        else:
            m = lists_cores[i-1] #take output of previous iteration as input
        patts, cores = replacecore_detect_homologue_cores(m, ru)
        lists_patts.append(patts)
        lists_cores.append(cores)
    #if there are any empty cores after all frag steps, get their index. Else, empty_cores_idx is itself empty list.
    empty_cores_idx = [i for i, j in enumerate(lists_cores[-1]) if j.GetNumAtoms() == 0]
    return lists_patts, lists_cores, empty_cores_idx

def my_ReplaceCore(mol,sidechain):
    '''Function that returns the original input mol if there is no matching sidechain present.'''
    rc = Chem.ReplaceCore(mol,sidechain)
    if rc is None: #mol does not contain sidechain whatsoever
        return mol
    else:
        return rc

def replacecore_longest_RU_match(mols, mat, ru):
    '''Function to replace the longest RU substructure match from each molecule with dummy atoms using RDKit's Chemical Transformations functionality. Returns remaining cores and the RU replaced by a dummy atom.'''
    cores = list()
    patts = list()
    for x,y in enumerate(mols):
        patts.append(ru[int(np.sum(mat[x])-1)])
        cores.append(my_ReplaceCore(y, patts[x]))
    return patts, cores

def detect_mols_made_of_ru(mols_with_ru, labels_mols_with_ru, empty_cores_idx):
    '''Function to detect, depict, then discard molecules made solely of RUs such as PEGs. Depictions written as png outputs.'''
    mols_made_of_ru = [j for i,j in enumerate(mols_with_ru) if (i in empty_cores_idx)] #isolate Mol and Label with empty core after first chopping i.e. entire mol is made of ru
    labels_made_of_ru = [j for i, j in enumerate(labels_mols_with_ru) if (i in empty_cores_idx)]
    mols_to_classify = [j for i,j in enumerate(mols_with_ru) if (i not in empty_cores_idx)]
    labels_to_classify = [j for i,j in enumerate(labels_mols_with_ru) if (i not in empty_cores_idx)]
    if len(mols_made_of_ru) > 0:
        pure_repeating_units = DrawMolsZoomed(mols_made_of_ru, labels_made_of_ru)
        pure_repeating_units.save("output_rmdum_tmf/mols_pure_repeating_units.png")
        print(str(len(mols_made_of_ru)) + " molecule(s) are made purely of repeating units of minimum length specified (default=3).")
    return mols_made_of_ru, labels_made_of_ru, mols_to_classify, labels_to_classify #includes non-series containing RU

def process_patts_cores(lists_patts, lists_cores, empty_cores_idx): #prepare dataframe???
    '''Filter out empty cores from patterns and cores lists in preparation for df assembly.'''
    
def generate_df_series(lists_patts, lists_cores etc.):

def DrawMolsZoomed(mols, legends, molsPerRow=3, subImgSize=(300, 300)):#, leg):
    """Function to draw rows of zoomed molecules. Credit Rocco Moretti."""
    nRows = len(mols) // molsPerRow
    if len(mols) % molsPerRow: nRows += 1
    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
    full_image = Image.new('RGB', fullSize )
    for ii, mol in enumerate(mols):
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)
        le = legends[ii]
        column = ii % molsPerRow
        row = ii // molsPerRow
        offset = ( column*subImgSize[0], row * subImgSize[1] )
        d2d = rdMolDraw2D.MolDraw2DCairo(subImgSize[0], subImgSize[1])
        d2d.DrawMolecule(mol,legend=le)
        d2d.FinishDrawing()
        sub = Image.open(BytesIO(d2d.GetDrawingText()))
        full_image.paste(sub,box=offset)
    return full_image

#from https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/CAHGTkV8sdfb4Q7FLn9C5MTwrqiJjHtnXK%2Bmz2SY3_4j2eAtevQ%40mail.gmail.com/#msg36477772

def delete_longest_RU_match(mols, mat, ru):
    '''Function to delete the longest RU substructure match from each molecule. Returns remaining cores and the RU deleted.'''
    cores = list()
    patt = list()
    for x,y in enumerate(mols):
        patt.append(ru[int(np.sum(mat[x])-1)])
        cores.append(AllChem.DeleteSubstructs(y, patt[x]))
    return patt, cores

def replace_longest_RU_match(mols, mat, ru):
    '''Function to replace the longest RU substructure match from each molecule with C. Returns remaining cores and the RU replaced by C.'''
    cores = list()
    patt = list()
    for x,y in enumerate(mols):
        patt.append(ru[int(np.sum(mat[x])-1)])
        cores.append(AllChem.ReplaceSubstructs(y, patt[x], Chem.MolFromSmiles('C'), replaceAll=True)[0])
    return patt, cores

def largest_core_molfrag_to_cano_smiles(cores2):
    '''Function to isolate the largest molecule fragment by NumAtoms remaining in the core object and convert to canonical SMILES.'''
    cores2_nonempty = [j for i, j in enumerate(cores2) if j.GetNumAtoms() > 0] # isolate non-empty cores after chopping
    #cores2_nonempty_molfrags = [Chem.GetMolFrags(co, asMols=True) for co in cores2_nonempty]  # cores2_nonempty_molfromsmarts
    #cores2_nonempty_largest_molfrag = []
    #for i, j in enumerate(cores2_nonempty_molfrags):  # get the largest molfrag of all molfrags
    #    cores2_nonempty_largest_molfrag.append(max(cores2_nonempty_molfrags[i],
    #                                               default=cores2_nonempty_molfrags[i],
    #                                               key=lambda m: m.GetNumAtoms()
    #                                               )
    #                                           )
    #cores2_nonempty_largest_molfrag_cano_smiles = [Chem.MolToSmiles(co, canonical=True) for co in cores2_nonempty_largest_molfrag]
    #cores2_nonempty_largest_molfrag_cano_smiles = [Chem.CanonSmiles(smi) for smi in cores2_nonempty_largest_molfrag_smiles]
    return cores2_nonempty#, cores2_nonempty_largest_molfrag_cano_smiles, cores2_nonempty_largest_molfrag

def GetNumFrags(mol):
    return(len(Chem.GetMolFrags(mol,asMols=True)))

def generate_output_summary(smiles_in, mols_classified, num_series, ru_in, mols_no_ru_matches, onememseries, mols_made_of_ru):
    mols_no_ru_matches = len(mols_no_ru_matches)
    onememseries = len(onememseries)
    mols_made_of_ru = len(mols_made_of_ru)
    with open("output_rmdum_tmf/classification-results.txt", "w") as text_file:
        text_file.write("Homologue classification for %s complete! \nRepeating Unit (RU) = %s.\n%s molecules were classified into %s homologous series.\n" % (smiles_in, ru_in, mols_classified, num_series))
        text_file.write("Molecules with no RU chains of min. three units in length: %s.\n" % (mols_no_ru_matches))
        text_file.write("Molecules having RU chains (min. three units) but are not members of any series: %s.\n" % (onememseries))
        text_file.write("Molecules consisting of purely RUs: %s." % (mols_made_of_ru))
