This repository provides tools to perform docking simulations with Autodock Vina for basic docking in a docking container. That means docking one ligand to a rigid receptor.

## Docking Container
We use AutoDock Vina software - O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking, with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010), 455-461, DOI 10.1002/jcc.21334.

The code for cantainer was prepared by Kamila Riedlová with the support of ChatGPT-4. After using this service, the author reviewed and edited the content as needed and take full responsibility for the content.

## Data for docking
Please upload the files (receptor and ligand) for which you wish to run the docking to the folder on your computer where you have the downloaded Dockerfile. In the Dockerfile, only modify the receptor and ligand names on these lines: 

    COPY receptor.pdb /data/
    
    COPY ligand.smi /data/

Furthermore, adjust the .json file for the receptor name, ligand, and also input the center and size parameters of the search space. 
    The center and size parameters define the search space for docking. Center_x,y,z parameters specify the x, y, and z coordinates of the search space's center. Ideally, this should be as close as possible to the expected interaction site between the ligand and receptor. If you are uncertain, you can open the receptor_file in the VMD program and find suitable pocket coordinates manually. For example: "x": 40.0, "y": 30.0, "z": 7.0.
    Size x, y, z determine the width, height, and depth of the search space. The search space should cover all possible interactions between the ligand and receptor but should also be small enough to minimize computational cost. Search space dimensions should be positive. For example: "x": 40.0, "y": 40.0, "z": 40.0.

## Build the container in the docking folder with the following commands: 
docker build -t name_of_image:1.0 .

    # For example: docker build -t docking:1.0 .
    
docker run -it --name name_of_container -v /home/username/your_folder/docking_folder:/data name_of_image:1.0 /bin/bash

    # For example: docker run -it --name docking-container -v /home/kiarka7/projects/docking:/data docking:1.0 /bin/bash

### Ligand Preparation
Acceptable formats for the ligand are SMILES (.smi), .mol2, .sdf, and .pdb. The script will automatically convert the .smi, .mol2, .sdf format to .pdb. If your ligand is already in the .pdb format, then use that. Subsequently, the .pdb file is converted to .pdbqt using Autodock Tools' prepare_ligand4.py. 

### Receptor Preparation
Acceptable format for the receptor is .pdb. The receptor is prepared using Autodock Tools' prepare_receptor4.py. Currently, it is possible to upload a receptor that contains atoms other than amino acids, especially crystal waters, ions, and compounds assisting crystallization. 

### Docking with AutoDock Vina
Once the .pdbqt formats for both the ligand and receptor are automatically generated, docking will automatically commence on these structures. AutoDock Vina software is utilized for the docking, ref.: O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking, with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010), 455-461, DOI 10.1002/jcc.21334.

### Output
In the folder where you run the container and where the input data resides, an 'output' folder will be generated. In it, you'll find a .zip file containing files for both the receptor and ligand in .pdbqt format, as well as the docking result as output.pdbqt.

### Analysis
Using the Vina forcefield, you should obtain an output.pdbqt with the best ligand docking score in term of negative binding affinity (kcal/mol). The first ligand structure in the output.pdbqt file should correspond to the best score. The output.pdbqt file can be visualized using Pymol tool (for both - output.pdbqt and receptor.pdbqt together), or you can also open it using AutoDockTools.
