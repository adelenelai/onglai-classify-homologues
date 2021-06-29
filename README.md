# A Cheminformatics Algorithm to Classify Homologous Series

## Introduction
Homologous series describe groups of chemical compounds containing a common core structure and a chain of growing repeating units.
 
 
 
 
 ## Requirements
 The homologues classification algorithm requires RDKit to be installed and the RDKit environment to be activated. This is best done using conda in the command line:
 
 ```
 conda create -c conda-forge -n my-rdkit-env rdkit
 
 conda activate my-rdkit-env
 ```
 
 More info here: https://www.rdkit.org/docs/Install.html




## Usage

```git clone https://github.com/adelenelai/classify_homologues
cd classify_homologues
conda activate my-rdkit-env       ## if not already activated
```

Then
```python classify_homologues.py <SMILES list> <Labels list> <repeating unit> 2>log```

Currently the only repeating unit is 'C', representing repeating alkyl groups i.e. -CH2-.

Try:

```
python classify_homols.py input/test1_smiles_23.csv input/test1_labels_23.csv C 2>log
```

Successful classification will generate an output directory containing .png plots of molecules, 1 plot per 1 classified series, and a .csv file containing machine-readable information on the classified series. 
Further optional outputs may be created depending on what the algorithm detects. Please see the sample output directories in this repository.



## License

This project is licensed under the Apache 2.0 License - see the [LICENSE](https://github.com/adelenelai/classify_homologues/blob/main/LICENSE) file for details.




## Authors

- [Adelene Lai](https://github.com/adelenelai)




## Acknowledgements
- 
## References

