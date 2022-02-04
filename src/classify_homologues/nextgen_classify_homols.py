import os
import sys
import rdkit.Chem as Chem
from utils import *
from rdkit.Chem import AllChem
from pathlib import Path
import rdkit.Chem.rdchem as rdchem
import rdkit.Chem.Draw as Draw
import rdkit.Chem.Descriptors
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole, MolsToGridImage
import pandas as pd
import numpy as np
import datamol as dm


print("Homologue classification started...")
#read in smiles and labels, prepare repeating units
smiles, mols = read_smiles_csv(sys.argv[1])
labels = read_labels_csv(sys.argv[2])
df = pd.DataFrame({ "SMILES":smiles, "Mols":mols, "Labels":labels})
ru_in = '[#6&H2]-'
ru = setup_repeating_unit(ru_in)
#tested '[#6&H2]-'
#tested '[#8]-[#6&H2]-[#6&H2]-'
#tested '[#6](-[#9])(-[#9])-'

#perform RU detection
mols_no_ru_matches, labels_mols_no_ru_matches, mols_with_ru, labels_mols_with_ru = detect_repeating_units(mols, labels, ru)

#perform core detection
patt1, cores1, patt2, cores2, empty_cores_idx = replacecore_detect_homologue_cores(mols_with_ru, ru)

#detect and output molecules made solely of RUs
Path("output").mkdir(parents=True, exist_ok=True)
mols_made_of_ru, labels_made_of_ru = detect_mols_made_of_ru(mols_with_ru, labels_mols_with_ru, empty_cores_idx)

#generate canonical SMILES of largest molecule fragment in cores
cores2_nonempty, cores2_nonempty_largest_molfrag_cano_smiles, cores2_nonempty_largest_molfrag = largest_core_molfrag_to_cano_smiles(cores2)


##construct dataframe for inspection
#finalise lists after filtering out mols_made_of_ru
mols_with_ru = [j for i,j in enumerate(mols_with_ru) if (i not in empty_cores_idx)]
labels_mols_with_ru = [j for i,j in enumerate(labels_mols_with_ru) if (i not in empty_cores_idx)]
#filter out row with empty core after first chopping from all cols and output
patt1 = [q for p,q in enumerate(patt1) if (p not in empty_cores_idx)]
cores1 = [q for p,q in enumerate(cores1) if (p not in empty_cores_idx)]
patt2 = [q for p,q in enumerate(patt2) if (p not in empty_cores_idx)]
cores2 = [q for p,q in enumerate(cores2) if (p not in empty_cores_idx)]




#build df to summarise SS match and chopping steps
imp_df_nonempty= pd.DataFrame({#"Smiles":smiles,
                   #"Mols":mols,
                    "Mols": mols_with_ru,
                    "Labels": labels_mols_with_ru,
                    "patt1": [j for i,j in enumerate(patt1) if cores2[i].GetNumAtoms()>0],
                    "Cores": [j for i,j in enumerate(cores1) if cores2[i].GetNumAtoms()>0],
                    "Cores2": cores2_nonempty,
                    "LargestMolFrag_sanitised": cores2_nonempty_largest_molfrag
                    })#"cano_smiles_cores2": cores2_cano_smiles



#populate new column in df with corresponding  cores2_nonempty_largestmolfrag_canosmiles for each row
imp_df_nonempty["canoSMILES_LargestMolfrag_sanitised"] = cores2_nonempty_largest_molfrag_cano_smiles

#assign SeriesNo to only those series with >1 member
imp_df_nonempty['SeriesNo'] = imp_df_nonempty.groupby('canoSMILES_LargestMolfrag_sanitised').filter(lambda group: len(group) > 1).groupby('canoSMILES_LargestMolfrag_sanitised').ngroup()


#merge imp_df_nonempty with df on mol
result = pd.merge(df, imp_df_nonempty, how="left", on=["Mols","Labels"])

#annotate no_alkyl mols AND series with 1 member to have negative SeriesNo
result['SeriesNo'] = result['SeriesNo'].fillna(-1)
result.SeriesNo = result.SeriesNo.astype(int)


#calc further identifiers and descr; write into log file
#with open('log.txt', 'w') as f:
 #   with redirect_stderr(f):
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

#depict final sanitised cores
if len(grpdmols.keys()) >0:
    final_cores = [Chem.MolFromSmiles(i) for i in grpdmols.keys()]
    leg_final_cores = [str(idx) for idx,y in enumerate(grpdmols.keys())]
    cores_summary = DrawMolsZoomed(final_cores, legends=leg_final_cores, molsPerRow=5)
    cores_summary.save("output/cores_summary.png")

grpdmols = [[Chem.MolFromSmiles(s) for s in g] for g in grpdmols]

list_grid_images = []
for i,j in enumerate(grpdmols):
    list_grid_images.append(DrawMolsZoomed(grpdmols[i], legends=lgs[i], molsPerRow=5))

#plot nans i.e. those with no alkyls and therefore not classified
#print(str(len(mols_no_alkyl_matches)) + ", " + str(len(labels_mols_no_alkyl_matches)))
if len(mols_no_ru_matches) > 0:
    nans = DrawMolsZoomed(mols=mols_no_ru_matches, legends=labels_mols_no_ru_matches, molsPerRow=5)
    nans.save("output/no_repeating_unit_matches.png")
#print(len(mols_no_alkyl_matches))
#nans = DrawMolsZoomed(mols=mols_no_alkyl_matches, legends=labels_mols_no_alkyl_matches, molsPerRow=5)
#nans.save("output/" + "no_alkyl_matches.png")

#save each plot per group
[img.save("output/" + str(idx) + ".png") for idx,img in enumerate(list_grid_images)]



#plot mols with alkyls but are 1-member series (i.e. not actually series)
#SeriesNo = -1, SMILES is not empty

onememseries = result.loc[(result['SeriesNo'] == -1) & (result['canoSMILES_LargestMolfrag_sanitised'].notnull())]
if len(onememseries.Mols) >0:
    mols_onememseries = [i for i in onememseries.Mols]
    labs_onememseries = [i for i in onememseries.Labels]
    pl_onememseries = DrawMolsZoomed(mols_onememseries,labs_onememseries,molsPerRow=5)
    pl_onememseries.save("output/non_series_containing_repeating_unit.png")
    print(str(len(onememseries.Mols))+ " molecule(s) have repeating unit matches of minimum x units but do not belong to any series.")


num_series = result.SeriesNo.max() + 1 #because zero-indexed
if num_series < 0:
    num_series= 0


mols_classified = len(result.Mols)-len(onememseries.Mols)-len(mols_no_ru_matches)-len(mols_made_of_ru)

print("Homologue classification complete! " + str(mols_classified) + " molecules have been classified into " +str(num_series) + " series." )

#output summary file
generate_output_summary(sys.argv[1], mols_classified, num_series, ru_in, mols_no_ru_matches, onememseries, mols_made_of_ru)
