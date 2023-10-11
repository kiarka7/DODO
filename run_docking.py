import os
import subprocess
import json
import argparse
import sys

MGLTOOLS_LIB_PATH = '/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs'
os.environ['PYTHONPATH'] = '/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs'

MGLTOOLS_BIN_PATH = '/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/'

def prepare_receptor(pdb_file, pdbqt_file):
    cmd = ["/usr/local/bin/python2.7", os.path.join(MGLTOOLS_BIN_PATH, 'prepare_receptor4.py'), '-r', pdb_file, '-o', pdbqt_file]
    subprocess.run(cmd)

def prepare_ligand(input_file, output_file):
    cmd = ["/usr/local/bin/python2.7", os.path.join(MGLTOOLS_BIN_PATH, 'prepare_ligand4.py'), '-l', input_file, '-o', output_file]
    current_directory = os.getcwd()
    os.chdir('/data')
    try:
        subprocess.run(cmd, check=True)
    finally:
        os.chdir(current_directory)

def convert_to_pdb(input_file, output_file):
    _, ext = os.path.splitext(input_file)
    if ext.lower() == '.smi':
        input_format = 'smi'
    elif ext.lower() == '.mol2':
        input_format = 'mol2'
    elif ext.lower() == '.cif':
        input_format = 'cif'
    elif ext.lower() == '.sdf':
        input_format = 'sdf'
    else:
        raise ValueError(f"Unsupported input format for conversion: {ext}")
    cmd = ["obabel", "-i", input_format, input_file, "-o", "pdb", "-O", output_file, "--gen3d"]
    subprocess.run(cmd, check=True)

DATA_PATH = "/data/"

def read_json(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    receptor = os.path.join(DATA_PATH, data["receptor"])
    ligand = os.path.join(DATA_PATH, data["ligand"])
    output = os.path.join(DATA_PATH, data["output"])
    
    center_x = data["center"]["x"]
    center_y = data["center"]["y"]
    center_z = data["center"]["z"]
    
    size_x = data["size"]["x"]
    size_y = data["size"]["y"]
    size_z = data["size"]["z"]
    
    return receptor, ligand, output, center_x, center_y, center_z, size_x, size_y, size_z

os.makedirs("/data/output", exist_ok=True)
os.chmod("/data/output", 0o777)

def check_file_exists(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found!")
    return filename

import zipfile

def create_zip(output_folder, files, zip_filename):
    with zipfile.ZipFile(zip_filename, 'w') as zipf:
        for file in files:
            zipf.write(file, os.path.basename(file))
    print(f"Saving to: {zip_filename}")


def main(json_file):
    receptor, ligand, output, center_x, center_y, center_z, size_x, size_y, size_z = read_json(json_file)

    receptor = check_file_exists(receptor)
    ligand = check_file_exists(ligand)

    # Convert receptor if it's not in pdbqt format
    base, ext = os.path.splitext(receptor)
    if ext.lower() == '.pdb':
        receptor_pdbqt = os.path.abspath(base + '.pdbqt')
        prepare_receptor(receptor, receptor_pdbqt)
        check_file_exists(receptor_pdbqt)
        receptor = receptor_pdbqt

    # Convert ligand if it's not in pdbqt format
    base, ext = os.path.splitext(ligand)
    if ext.lower() == '.pdb':
        ligand_pdbqt = os.path.abspath(base + '.pdbqt')
        print(f"Converting {ligand} to {ligand_pdbqt} ...")
        prepare_ligand(ligand, ligand_pdbqt)
        check_file_exists(ligand_pdbqt)
        ligand = ligand_pdbqt
    elif ext.lower() in ['.smi', '.mol2', '.cif', '.sdf']:
        ligand_pdb = os.path.join('/data', base + '.pdb')
        print(f"Converting {ligand} to {ligand_pdb} ...")
        convert_to_pdb(ligand, ligand_pdb)
        check_file_exists(ligand_pdb)
        ligand_pdbqt = os.path.abspath(base + '.pdbqt')
        print(f"Converting {ligand_pdb} to {ligand_pdbqt} ...")
        prepare_ligand(ligand_pdb, ligand_pdbqt)
        check_file_exists(ligand_pdbqt)
        ligand = ligand_pdbqt
    elif ext.lower() != '.pdbqt':
        raise ValueError("Unsupported ligand file format: " + ext)
    
    # Run AutoDock Vina
    vina_cmd = [
        'vina',
        '--receptor', receptor,
        '--ligand', ligand,
        '--out', output,
        '--center_x', str(center_x),
        '--center_y', str(center_y),
        '--center_z', str(center_z),
        '--size_x', str(size_x),
        '--size_y', str(size_y),
        '--size_z', str(size_z),
    ]
    subprocess.run(vina_cmd)

    output_folder = os.path.join(DATA_PATH, "output")
    os.makedirs(output_folder, exist_ok=True)
    zip_filename = os.path.join(output_folder, "results.zip")
    files_to_zip = [output, receptor, ligand]  # přidán ligand do zipu
    create_zip(output_folder, files_to_zip, zip_filename)

    print(f"Results saved to: {zip_filename}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 docking_script.py path_to_json_file")
        sys.exit(1)
    json_file = sys.argv[1]
    main(json_file)