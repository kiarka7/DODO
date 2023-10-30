# DODO - Docking in docking container 
This repository provides tools to perform docking simulations with [Autodock Vina](https://vina.scripps.edu/) for basic docking in a docking container. That means docking one ligand to a rigid receptor.

## Usage

In order to carry out the docking, you will need a receptor and a ligand file (see below for valid file formats) and a json file with parameters. All these three files need to reside in the same directory (the receptor and ligand file can actually be in a subdirectory of the directory where the params file resides). This directory is then mapped to a `data` directory in the Docker image and after the docking, the output zip file is stored there.


### Params file

The name of the params file needs to be `docking_parameters.json`. It is a json file, that, besides the input files, also defines the center and size parameters of the search space determining the search space for docking. 

The _x_, _y_, and _z_ coordinates of the search space's center  should be as close as possible to the expected interaction site between the ligand and receptor. If you are uncertain, you can open the receptor_file in the [VMD](https://www.ks.uiuc.edu/Research/vmd/) or [PyMOL](https://pymol.org/2/) and find suitable pocket coordinates manually. The Size _x_, _y_, _z_ parameters determine the width, height, and depth of the search space. The search space should cover all possible interactions between the ligand and receptor but should also be small enough to minimize computational cost. Search space dimensions should be positive. 

See an example of a params file in the test_files folder, eg 1za1_D.

### Running the continer 

To run the contain, you first need to build the image.

```console
docker build -t docking .
```

Once the image is ready, you can create a container from the built image. Follows an example with the one of the examples from the [test_files](test_files) folder.

```console
 docker run -it --name docking-container -v /mnt/c/projects/git/docking-docker/test_files/1za1_D/:/data docking
```

The outputs will be stored in a zip file in the mapped directory.
   
## Docking process

### Ligand Preparation
Acceptable formats for the ligand are SMILES `.smi`, `.mol2`, `.sdf`, and `.pdb`. The script will automatically convert the `.smi`, `.mol2`, `.sdf` format to `.pdb`. If your ligand is already in the `.pdb` format, the script uses that one directly. Subsequently, the `.pdb` file is converted to `.pdbqt` using Autodock Tools' [`prepare_ligand4.py`](https://github.com/sahrendt0/Scripts/blob/master/docking_scripts/prepare_ligand4.py). 

### Receptor Preparation
Acceptable format for the receptor is `.pdb`. The receptor is prepared using Autodock Tools' [`prepare_receptor4.py`](https://github.com/sahrendt0/Scripts/blob/master/docking_scripts/prepare_receptor4.py). Currently, it is possible to upload a receptor that contains atoms other than amino acids, especially crystal waters, ions, and compounds assisting crystallization. 

### Docking with AutoDock Vina
Once the `.pdbqt`` formats for both the ligand and receptor are automatically generated, docking will automatically commence on these structures using the  AutoDock Vina software.

### Output
In the folder where you run the container and where the input data resides, an 'output' folder will be generated. In it, you'll find a `.zip` file containing files for both the receptor and ligand in `.pdbqt` format, as well as the docking result as `output.pdbqt`.

### Analysis
Using the Vina forcefield, you should obtain an `output.pdbqt` with the best ligand docking score in term of negative binding affinity (kcal/mol). The first ligand structure in the `output.pdbqt` file should correspond to the best score. The output.pdbqt file can be visualized using Pymol tool (for both - `output.pdbqt` and `receptor.pdbqt` together), or you can also open it using AutoDockTools.

## Acknowledgements

- O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking, with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010), 455-461, [DOI: 10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334).
