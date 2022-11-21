# An Algorithm to Classify Homologous Series
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/adelenelai/classify_homologues/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/adelenelai/onglai-classify-homologues.svg)](https://GitHub.com/adelenelai/onglai-classify-homologues/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/adelenelai/onglai-classify-homologues.svg)](https://GitHub.com/adelenelai/onglai-classify-homologues/graphs/contributors/)
[![DOI](https://zenodo.org/badge/381339802.svg)](https://zenodo.org/badge/latestdoi/381339802)
[![GitHub release](https://img.shields.io/github/release/adelenelai/onglai-classify-homologues.svg)](https://github.com/adelenelai/onglai-classify-homologues/releases/)
[![PyPI version fury.io](https://badge.fury.io/py/onglai.svg)](https://pypi.python.org/pypi/onglai/)



## Introduction
Homologous series are groups of chemical compounds sharing the same core structure(s) and different numbers of repeating units (RU) connected end-to-end.

This is an open-source algorithm to classify homologous series within compound datasets provided as SMILES, implemented using the RDKit.

For example, these series were classified in [COCONUT](https://coconut.naturalproducts.net/) and the [NORMAN Suspect List Exchange](https://www.norman-network.com/nds/SLE/), datasets containing natural products and environmental chemicals respectively.


CH2 Repeating Unit:
![coconut-hs](https://github.com/adelenelai/onglai-classify-homologues/blob/main/5027.png)

CF2 Repeating Unit:
![norman-hs](https://github.com/adelenelai/onglai-classify-homologues/blob/main/11_epoxy.png)



## Requirements
 The algorithm requires RDKit to be [installed](https://www.rdkit.org/docs/Install.html) via `conda-forge`.

 ```shell
 $ conda create -c conda-forge -n my-rdkit-env rdkit
 $ conda activate my-rdkit-env
 ```


## Installation

```shell
$ git clone https://github.com/adelenelai/onglai-classify-homologues
$ cd classify_homologues
$ pip install -e .
```
Note that pip installing the package is not enough; in addition, the repo must be cloned from GitHub because the algorithm runs as a script (see below).

Alternatively:
```python
#from PyPI
$ pip install onglai-classify-homologues
```

## Usage

Run:

```shell
$ python nextgen_classify_homols.py [-in <arg>] [-s <arg>] [-n <arg>] [-ru <arg>] [-min <arg>] [-max <arg>] 2>log
```

| Flag | Description |
| --- | ----------- |
| -in --input_csv <arg> | path to input CSV containing SMILES and Name columns|
| -s --smiles <arg> | name of column containing SMILES. Default is SMILES.|
| -n --names <arg> | name of column containing Names. Default is Name.|
| -ru --repeatingunits <arg> | chemical RU as SMARTS, enclosed within speech marks. Default is CH2 i.e., '[#6&H2]'. |
| -min --min_RU_in <arg> | minimum length of RU chain, default is 3|
| -max --max__RU_in <arg> | maximum length of RU chain, default is 30 |
| -f --frag_steps <arg> | no. times to fragment molecules to obtain cores, the default is 2 |


Try:
```shell
$ python nextgen_classify_homols.py -in ../../tests/test1_23.csv -s SMILES -n Name -ru '[#6&H2]' -min 3 -max 30 -f 2 2>log
```

Successful classification will generate an `output` directory containing the following files:

1. A TXT file containing the summary of classification results and explanation of outputs (series_no codes)
2. A CSV file containing 8 columns: `series_no`, `cpd_name`, `CanoSmiles_FinalCores`, `SMILES`, `InChI`, `InChIKey`, `molecular_formula` and `monoisotopic_mass`. The first column `series_no` contains the results of the homologous series classification. `CanoSmiles_FinalCores` indicates the common core shared by all members within a given series. The remaining columns contain information calculated based on the `SMILES`.
3. A TXT file of unparseable SMILES that were removed (if all SMILES were parsed OK, then empty)


### Reproducing Classification described in Lai et al.

Classification using default settings as described above. Code below runs for sample datasets provided in `input/`, full datasets have been archived on [Zenodo](https://doi.org/10.5281/zenodo.6958826) (amend `-in` accordingly to classify full datasets).

```
#activate your rdkit environment

#NORMAN-SLE
$ python nextgen_classify_homols.py -in ../../input/pubchem_norman_sle_tree_parentcid_98116_2022-03-21_from115115_trial.csv -s isosmiles -n cmpdname 2>log

#PubChemLite
$ python nextgen_classify_homols.py -in ../../input/PubChemLite_exposomics_20220225_trial.csv -n CompoundName 2>log

#COCONUT
$ python nextgen_classify_homols.py -in ../../input/COCONUT_DB_2021-11_trial.txt 2>log
```



## References and Links
* Lai, A., Schaub, J., Steinbeck, C., Schymanski, E. L. An Algorithm to Classify Homologous Series in Compound Datasets. [Preprint](https://doi.org/10.21203/rs.3.rs-2019306/v1)
* [Poster](https://zenodo.org/record/6491204) presented at the 17th German Cheminformatics Conference, Garmisch-Partenkirchen, Germany (May 8-10, 2022)


## License

This project is licensed under Apache 2.0  - see [LICENSE](https://github.com/adelenelai/onglai-classify-homologues/blob/main/LICENSE) for details.


## Our Research Groups
[Environmental Cheminformatics Group](https://wwwen.uni.lu/lcsb/research/environmental_cheminformatics) at the


[<p align="center"><img src="https://github.com/adelenelai/onglai-classify-homologues/blob/main/logo_LCSB_UL.png" width='50%'></p>](https://wwwen.uni.lu/lcsb)

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
