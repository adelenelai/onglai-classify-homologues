# A Cheminformatics Algorithm to Classify Homologous Chemical Series

## Introduction
Homologous chemical series describe groups of chemical compounds containing common core structure(s) and a chain of growing repeating units (RUs). This is an open-source algorithm to classify homologous series implemented in the RDKit.



## Requirements
 The homologues classification algorithm requires RDKit to be [installed](https://www.rdkit.org/docs/Install.html) and the RDKit environment to be activated. This is best done using conda in the command line:

 ```
 $ conda create -c conda-forge -n my-rdkit-env rdkit
 $ conda activate my-rdkit-env
 $ git clone https://github.com/adelenelai/classify_homologues
 $ cd classify_homologues
 ```



## Usage
The algorithm runs in the command line interface as below:

```
python classify_homologues.py [-s <arg>] [-l <arg>] [-ru <arg>] [-min <arg>] [-max <arg>] 2>log
```

| Flag | Description |
| --- | ----------- |
| -s --smiles <arg> | path to CSV list of SMILES |
| -l --labels <arg> | path to CSV list of Labels (molecule names) |
| -ru --repeatingunits <arg> | chemical RU as SMARTS, enclosed within speech marks. Default is CH2 i.e., '[#6&H2]'. |
| -min --min_RU_in <arg> | minimum length of RU chain, default is 3|
| -max --max__RU_in <arg> | maximum length of RU chain, default is 30 |
| -f --frag_steps <arg> | no. times to fragment molecules to obtain cores, default is 2 |


Try:
```
python nextgen_classify_homols.py -s input/test1_smiles_23.csv -l input/test1_labels_23.csv -ru '[#6&H2]' -min 3 -max 5 -f 3 2>log
```

Successful classification will generate an output directory containing .png plots of molecules (e.g., 1  plot per 1 classified series), a CSV file containing machine-readable information on the classified series, and a TXT file containing a classification summary.


Further optional outputs may be created depending on what the algorithm detects (e.g. molecules with no repeating units detected). Please see the sample output directories in this repository e.g., `output_test1/` for the example given above.



## License

This project is licensed under the Apache 2.0 License - see the [LICENSE](https://github.com/adelenelai/classify_homologues/blob/main/LICENSE) file for details.




## Authors

- [Adelene Lai](https://github.com/adelenelai)



## Acknowledgements

-
## References

## Our Research Group
[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
