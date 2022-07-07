# An Algorithm to Classify Homologous Series

## Introduction


## Requirements
 The algorithm requires RDKit to be [installed](https://www.rdkit.org/docs/Install.html) and the RDKit environment to be activated. This is best done using conda in the command line:

 ```
 $ conda create -c conda-forge -n my-rdkit-env rdkit
 $ conda activate my-rdkit-env
 ```

## Installation
The code can be installed directly from GitHub with:

```shell
$ pip install git+https://github.com/adelenelai/classify_homologues.git
```


Or in editable mode with:

```shell
$ git clone https://github.com/adelenelai/classify_homologues
$ cd classify_homologues
$ pip install -e .
```

## Usage
```
python nextgen_classify_homols.py [-s <arg>] [-l <arg>] [-ru <arg>] [-min <arg>] [-max <arg>] 2>log
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

Successful classification will generate an `output` directory containing the following files:

1. a TXT file containing the summary of classification results
2. a CSV file containing 8 columns: `series_no`, `cpd_name`, `CanoSmiles_FinalCores`, `SMILES`, `InChI`, `InChIKey`, `molecular_formula` and `monoisotopic_mass`. The first column `series_no` contains the results of the homologous series classification. `CanoSmiles_FinalCores` indicates the common core shared by all members within a given series.  Columns `SMILES` and `cpd_name` were the original inputs to the `-s` and `-l` flags respectively. The remaining columns contain information calculated based on the `SMILES`.
3. a TXT file of unparseable SMILES that were removed (if all SMILES were parsed OK, then empty)


## References and Links
* *publication coming soon*
* [Poster](https://zenodo.org/record/6491204) presented at the 17th German Cheminforamtics Conference, Garmisch-Partenkirchen, Germany (May 8-10, 2022)


## License

This project is licensed under the Apache 2.0 License - see the [LICENSE](https://github.com/adelenelai/classify_homologues/blob/main/LICENSE) file for details.




## Authors

- [Adelene Lai](https://github.com/adelenelai)



## Acknowledgements

-

## Our Research Group
[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
