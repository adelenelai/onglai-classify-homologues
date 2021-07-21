import os
import sys
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import rdkit.Chem.rdchem as rdchem
import rdkit.Chem.Draw as Draw
import rdkit.Chem.Descriptors
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole, MolsToGridImage
PandasTools.RenderImagesInAllDataFrames(images=True)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit.Chem import rdFMCS
from itertools import compress
from collections import defaultdict, OrderedDict, Counter
from timeit import default_timer as timer
from contextlib import redirect_stderr

#from https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/CAHGTkV8sdfb4Q7FLn9C5MTwrqiJjHtnXK%2Bmz2SY3_4j2eAtevQ%40mail.gmail.com/#msg36477772
from rdkit.Chem.Draw import rdMolDraw2D
try:
    import Image
except ImportError:
    from PIL import Image
from io import BytesIO

def DrawMolsZoomed(mols, legends, molsPerRow=3, subImgSize=(300, 300)):#, leg): #https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDraw2D
    nRows = len(mols) // molsPerRow
    if len(mols) % molsPerRow: nRows += 1
    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
    full_image = Image.new('RGBA', fullSize )
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


Chem.WrapLogs()
#from io import StringIO
#sio = sys.stderr = StringIO()


#mols
with open(sys.argv[1]) as f:
    smiles = [line.strip().replace('"','') for line in f]

mols = [AllChem.MolFromSmiles(smile) for smile in smiles]

#add explicitHs to terminal atoms only (for dummy atom matching)
mols = [Chem.AddHs(m,onlyOnAtoms=(m.GetNumAtoms()-1,0)) for m in mols]
mols = [Chem.AddHs(m,onlyOnAtoms=(0,0)) for m in mols]

#labels
with open(sys.argv[2]) as f:
    labels = [line.strip().replace('"','') for line in f]

df = pd.DataFrame({ "SMILES":smiles,
                     "Mols":mols,
                    "Labels":labels
                  })

#read in alkyl chains
#with open('alkylchains.txt') as p:
#    smiles_alkylchains = [line.strip() for line in p]

#repeating unit
smiles_ru = []
#enumerate
[smiles_ru.append(sys.argv[3]*i) for i in range(1,31)]




#see add_dummy_atoms_bookend.ipynb as of 2021-07-08
#add start and end dummy atoms
smiles_ru = ['*' + ru + '*' for ru in smiles_ru]

ru = [AllChem.MolFromSmiles(smi) for smi in smiles_ru] #includes C and CC

minlen = int(sys.argv[4])-1 if len(sys.argv) >= 4 else 2

#prepare to generate MergeQueryHs
ru_explicitH = [Chem.AddHs(r) for r in ru[minlen:]]
#ru_explicitH = [Chem.AddHs(r,onlyOnAtoms=(range(1,r.GetNumAtoms()-1))) for r in ru[2:]] #automatically discards meth- and eth-


#with implicitHs i.e. MergeQueryHs
ru_implicit_queryH = [Chem.MergeQueryHs(m) for m in ru_explicitH]


#make dummy atoms queries i.e. [0] to *
ru_implicit_queryH = [Chem.AdjustQueryProperties(ru) for ru in ru_implicit_queryH]

#initialise foo - the 2d array showing whether there are alkyl substructure matches (1 or 0)
foo = np.zeros((len(mols), len(ru_implicit_queryH)))

for a,b in enumerate(mols):
    for i,j in enumerate(ru_implicit_queryH):
        foo[a,i] = mols[a].HasSubstructMatch(j) #gives 1s and 0s in alkyl chain matching

#analyse foo_array_sums
foo_array_sums = []
length_ru_chains_to_chop = []
for x,y in enumerate(foo):
    foo_array_sums.append(int(np.sum(foo[x])))
    length_ru_chains_to_chop.append(foo_array_sums[x] + 2) # n = sum + 2, where n is the length of C chain to chop

#how many have no alkyl chains? i.e. value = 0
n_mols_no_ru = foo_array_sums.count(0)
#print(str(n_mols_no_ru) + " mols have no repeating unit chains of minimum " + str(sys.argv[4]) + " repeating units in length.")
#remove mols with no alkyl matches
fil_ru = []
fil_ru = [bool(x) for x in foo_array_sums] #those which are False have array_sum = 0 i.e. no alkyls
mols_with_ru = list(compress(mols,fil_ru))
labels_mols_with_ru = list(compress(labels,fil_ru))
#len(mols_with_alkyls)

mols_no_ru_matches = list(compress(mols, [not i for i in fil_ru]))
mols_no_ru_matches = [Chem.RemoveHs(m) for m in mols_no_ru_matches] #for plotting
labels_mols_no_ru_matches = list(compress(labels, [not i for i in fil_ru]))
#len(mols_no_alkyl_matches)

#redo chopping and dictionary generation using mols list: mols_with_alkyls
#initialise 2d array to do longest alkyl chain matching
foo = np.zeros((len(mols_with_ru), len(ru_implicit_queryH)))
for a,b in enumerate(mols_with_ru):
    for i,j in enumerate(ru_implicit_queryH):
        foo[a,i] = mols_with_ru[a].HasSubstructMatch(j) #1s and 0s in alkyl chain matching

#start 2-step removal in one go without generating new idxs
##
#Step 1
cores = list()
patt1 = list()
for x,y in enumerate(mols_with_ru):
    patt1.append(ru_implicit_queryH[int(np.sum(foo[x])-1)]) #cannot use extend; Mol object is not iterable
    cores.append(AllChem.DeleteSubstructs(y,patt1[x]))

#initialise 2d array for second removal
goo = np.zeros((len(mols_with_ru), len(ru_implicit_queryH)))

for a,b in enumerate(cores):
    for i,j in enumerate(ru_implicit_queryH):
        goo[a,i] = cores[a].HasSubstructMatch(j)

#2-step removal in one go without generating new idxs
##
#Step 2
cores2 = list()
patt2 = list()
for x,y in enumerate(cores):
    patt2.append(ru_implicit_queryH[int(np.sum(goo[x])-1)]) #if np.sum is 0, could have issues
    cores2.append(AllChem.DeleteSubstructs(y,patt2[x]))

#isolate empty cores after first chopping
#happens when mol is 100% made of RU
empty_cores_idx = [i for i,j in enumerate(cores) if j.GetNumAtoms()==0]

#isolate Mol and Label with empty core after first chopping i.e. entire mol is made of ru
mols_made_of_ru = [j for i,j in enumerate(mols_with_ru) if (i in empty_cores_idx)]
mols_made_of_ru = [Chem.RemoveHs(m) for m in mols_made_of_ru] #remove explicitHs for plotting
labels_made_of_ru = [j for i,j in enumerate(labels_mols_with_ru) if (i in empty_cores_idx)]

#finalise lists after filtering out mols_made_of_ru
mols_with_ru = [j for i,j in enumerate(mols_with_ru) if (i not in empty_cores_idx)]
labels_mols_with_ru = [j for i,j in enumerate(labels_mols_with_ru) if (i not in empty_cores_idx)]

#isolate non-empty cores after chopping
cores2_nonempty = [j for i,j in enumerate(cores2) if j.GetNumAtoms()>0]

 #make output dir
os.makedirs("output"+"/")
if len(mols_made_of_ru) >0:
    pure_repeating_units = DrawMolsZoomed(mols_made_of_ru,labels_made_of_ru)
    pure_repeating_units.save("output/mols_pure_repeating_units.png")
    print(str(len(mols_made_of_ru)) + " molecule(s) are made purely of repeating units of minimum length " + str(sys.argv[4] + "."))

#filter out row with empty core after first chopping from all cols and output
patt1 = [q for p,q in enumerate(patt1) if (p not in empty_cores_idx)]
cores = [q for p,q in enumerate(cores) if (p not in empty_cores_idx)]
patt2 = [q for p,q in enumerate(patt2) if (p not in empty_cores_idx)]
cores2 = [q for p,q in enumerate(cores2) if (p not in empty_cores_idx)]



#first convert to smarts, then molfromsmarts, then sanitize
cores2_nonempty_smarts = [Chem.rdmolfiles.MolToSmarts(co) for co in cores2_nonempty]
cores2_nonempty_molfromsmarts = [Chem.MolFromSmarts(co) for co in cores2_nonempty_smarts] #generates QueryAtoms
[Chem.SanitizeMol(co) for co in cores2_nonempty_molfromsmarts]
#get molfrags
cores2_nonempty_molfrags = [Chem.GetMolFrags(co, asMols=True) for co in cores2_nonempty_molfromsmarts]
#get the largest molfrag of all molfrags
cores2_nonempty_largest_molfrag = []
for i,j in enumerate(cores2_nonempty_molfrags):
    cores2_nonempty_largest_molfrag.append(max(cores2_nonempty_molfrags[i], default=cores2_nonempty_molfrags[i], key=lambda m:m.GetNumAtoms()))

#build df to summarise SS match and chopping steps
imp_df_nonempty= pd.DataFrame({#"Smiles":smiles,
                   #"Mols":mols,
                    "Mols": mols_with_ru,
                    "Labels": labels_mols_with_ru,
                    "patt1": [j for i,j in enumerate(patt1) if cores2[i].GetNumAtoms()>0],
                    "Cores": [j for i,j in enumerate(cores) if cores2[i].GetNumAtoms()>0],
                    "patt2": [j for i,j in enumerate(patt2) if cores2[i].GetNumAtoms()>0],
                    "Cores2": cores2_nonempty,
                    "LargestMolFrag_sanitised": cores2_nonempty_largest_molfrag
                    })#"cano_smiles_cores2": cores2_cano_smiles

#get canonical smiles of largest molfrag
cores2_nonempty_largest_molfrag_smiles = [Chem.MolToSmiles(co) for co in cores2_nonempty_largest_molfrag]
cores2_nonempty_largest_molfrag_cano_smiles = [Chem.CanonSmiles(smi) for smi in cores2_nonempty_largest_molfrag_smiles]

#populate new column in df with corresponding  cores2_nonempty_largestmolfrag_canosmiles for each row
imp_df_nonempty["canoSMILES_LargestMolfrag_sanitised"] = cores2_nonempty_largest_molfrag_cano_smiles

#assign SeriesNo to only those series with >1 member
imp_df_nonempty['SeriesNo'] = imp_df_nonempty.groupby('canoSMILES_LargestMolfrag_sanitised').filter(lambda group: len(group) > 1).groupby('canoSMILES_LargestMolfrag_sanitised').ngroup()


#merge imp_df_nonempty with df on mol
result = pd.merge(df, imp_df_nonempty, how="left", on=["Mols","Labels"])

#annotate no_alkyl mols AND series with 1 member to have negative SeriesNo
result['SeriesNo'] = result['SeriesNo'].fillna(-1)
result.SeriesNo = result.SeriesNo.astype(int)


#calc further identifiers and descr
with open('yourfile.txt', 'w') as f:
    with redirect_stderr(f):
        inchis = [Chem.inchi.MolToInchi(i) for i in result.Mols]
        inchikeys = [Chem.inchi.MolToInchiKey(i) for i in result.Mols]
        mf = [Chem.rdMolDescriptors.CalcMolFormula(i) for i in result.Mols]
        monoiso_mass = [round(Chem.Descriptors.ExactMolWt(i),4) for i in result.Mols]

out = result[["SeriesNo","Labels","canoSMILES_LargestMolfrag_sanitised","SMILES", ]].copy()
out['InChI'] = inchis
out['InChIKey'] = inchikeys
out['molecular_formula'] = mf
out['monoisotopic_mass'] = monoiso_mass
out.rename(columns={"SeriesNo":"series_no", "Labels":"series_name", "canoSMILES_LargestMolfrag_sanitised": "common_core"},inplace=True)
out.to_csv('output/' + 'classified_series.csv',index=False)

#plots per group - only series which have >1 member i.e. actual series
result_pos_serno = result[result["SeriesNo"] > -1]
#legends
lgs = [i for i in result_pos_serno.groupby('canoSMILES_LargestMolfrag_sanitised').Labels.apply(list)]
for i,j in enumerate(lgs):
    lgs[i] = lgs[i] + ["core"]

grpdmols = result_pos_serno.groupby('canoSMILES_LargestMolfrag_sanitised').SMILES.apply(list) #lists of SMILES
for i,j in enumerate(grpdmols):
    grpdmols[i] = grpdmols[i] + [grpdmols.keys()[i]]

grpdmols = [[Chem.MolFromSmiles(s) for s in g] for g in grpdmols]

list_grid_images = []
for i,j in enumerate(grpdmols):
    list_grid_images.append(DrawMolsZoomed(grpdmols[i],legends=lgs[i],molsPerRow=5))

#plot nans i.e. those with no alkyls and therefore not classified
#print(str(len(mols_no_alkyl_matches)) + ", " + str(len(labels_mols_no_alkyl_matches)))
if len(mols_no_ru_matches) > 0:
    nans = DrawMolsZoomed(mols=mols_no_ru_matches, legends=labels_mols_no_ru_matches, molsPerRow=5)
    nans.save("output/" + "no_repeating_unit_matches.png")
    print(str(len(mols_no_ru_matches)) + " molecule(s) have no repeating unit matches of minimum " + str(sys.argv[4]) + " units.")
#nans = DrawMolsZoomed(mols=mols_no_alkyl_matches, legends=labels_mols_no_alkyl_matches, molsPerRow=5)
#nans.save("output/" + "no_alkyl_matches.png")

#save each plot per group
[img.save("output/" + str(idx) + ".png") for idx,img in enumerate(list_grid_images)]


#plot mols with alkyls but are 1-member series (i.e. not actually series)
#SeriesNo = -1, SMILES is not empty

onememseries = result.loc[(result['SeriesNo'] == -1) & (result['canoSMILES_LargestMolfrag_sanitised'].notnull())]
if len(onememseries.Mols) >0:
    mols_onememseries = [i for i in onememseries.Mols]
    mols_onememseries = [Chem.RemoveHs(m) for m in mols_onememseries]
    labs_onememseries = [i for i in onememseries.Labels]
    pl_onememseries = DrawMolsZoomed(mols_onememseries,labs_onememseries,molsPerRow=5)
    pl_onememseries.save("output/non_series_containing_repeating_unit.png")
    print(str(len(onememseries.Mols))+ " molecule(s) have repeating unit matches of minimum " + str(sys.argv[4]) + " units but do not belong to any series.")

num_series = result.SeriesNo.max()
#num_series = num_series+1

#if num_series < 0:
#    num_series= 0

#sys.stderr = old_stderr
#assert sio.getvalue() != ""
mols_classified = len(result.Mols)-len(onememseries.Mols)-len(mols_no_ru_matches)-len(mols_made_of_ru)
print("Homologue classification complete! " + str(mols_classified) + " molecules have been classified into "+ str(num_series+1) + " series." )
