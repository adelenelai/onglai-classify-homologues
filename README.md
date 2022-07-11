# An Algorithm to Classify Homologous Series
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/adelenelai/classify_homologues/graphs/commit-activity)
[![DOI](https://zenodo.org/badge/381339802.svg)](https://zenodo.org/badge/latestdoi/381339802)



## Introduction
Homologous series are groups of chemical compounds sharing the same core structure(s) and different numbers of repeating units (RU) connected end-to-end.

This is an open-source algorithm to classify homologous series within compound datasets provided as SMILES, implemented using the RDKit. 

For example, these series were classified in [COCONUT](https://coconut.naturalproducts.net/) and the [NORMAN Suspect List Exchange](https://www.norman-network.com/nds/SLE/).


CH2 Repeating Unit:
![GitHub Logo](https://github.com/adelenelai/classify_homologues/blob/main/5027.png)

CF2 Repeating Unit:
![GitHub Logo](https://github.com/adelenelai/classify_homologues/blob/main/11_epoxy.png)



## Requirements
 The algorithm requires RDKit to be [installed](https://www.rdkit.org/docs/Install.html) via `conda-forge`.

 ```
 $ conda create -c conda-forge -n my-rdkit-env rdkit
 $ conda activate my-rdkit-env
 ```


## Installation

```shell
$ git clone https://github.com/adelenelai/classify_homologues
$ cd classify_homologues
$ pip install -e .
```
Note that pip installing the package is not enough; in addition, the repo must be cloned from GitHub because the algorithm runs as a script (see below).

## Usage

Run:

```
$ python nextgen_classify_homols.py [-s <arg>] [-l <arg>] [-ru <arg>] [-min <arg>] [-max <arg>] 2>log
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
$ python src/classify_homologues/nextgen_classify_homols.py -s tests/test1_smiles_23.csv -l tests/test1_labels_23.csv -ru '[#6&H2]' -min 3 -max 5 -f 3 2>log
```

Successful classification will generate an `output` directory containing the following files:

1. a TXT file containing the summary of classification results
2. a CSV file containing 8 columns: `series_no`, `cpd_name`, `CanoSmiles_FinalCores`, `SMILES`, `InChI`, `InChIKey`, `molecular_formula` and `monoisotopic_mass`. The first column `series_no` contains the results of the homologous series classification. `CanoSmiles_FinalCores` indicates the common core shared by all members within a given series.  Columns `SMILES` and `cpd_name` were the original inputs to the `-s` and `-l` flags respectively. The remaining columns contain information calculated based on the `SMILES`.
3. a TXT file of unparseable SMILES that were removed (if all SMILES were parsed OK, then empty)


## References and Links
* *publication coming soon 2022!*
* [Poster](https://zenodo.org/record/6491204) presented at the 17th German Cheminforamtics Conference, Garmisch-Partenkirchen, Germany (May 8-10, 2022)


## License

This project is licensed under Apache 2.0  - see the [LICENSE](https://github.com/adelenelai/classify_homologues/blob/main/LICENSE) file for details.


## Our Research Groups
Environmental Cheminformatics Group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
