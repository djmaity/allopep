# allopep
AlloPep: A Pipeline for De Novo Design of Allosteric Peptides

## Installation
```
bash install.sh
```

Run the following command to create the conda environment with the name `allopep` using the `environment.yaml` file
```
conda env create -f environment.yml -n allopep
```

Ensure that the `allopep` conda environment is acitvated before
run the following command from the allopep directory

## Running
```
python allopep.py input/structure.pdb
```
## Output

The output will be placed in the output directory created within the allopep directory

Please make sure the input file, e.g., structure.pdb is in not placed in the output directory. APOP creates a structure.pdb file in the output directory after removing water and hetero atoms, which will overwirte the input structure.pdb if placed in the output directory.
