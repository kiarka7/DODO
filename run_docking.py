import os
import subprocess
import json
import sys
import zipfile
import shutil

MGLTOOLS_LIB_PATH = '/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs'
MGLTOOLS_BIN_PATH = '/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/'
DATA_PATH = "/data/"
WORK_PATH ="/tmp/"

def filter_pdb(input_pdb, output_pdb):
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                f_out.write(line)

    print(f"Filtered PDB saved to {output_pdb}")

def prepare_receptor(pdb_file, pdbqt_file):
    cmd = ["/usr/local/bin/python2.7", os.path.join(MGLTOOLS_BIN_PATH, 'prepare_receptor4.py'), '-r', pdb_file, '-o', pdbqt_file]
    subprocess.run(cmd, check=True)

def prepare_ligand(input_file, output_file):
    cmd = ["/usr/local/bin/python2.7", os.path.join(MGLTOOLS_BIN_PATH, 'prepare_ligand4.py'), '-l', input_file, '-o', output_file]
    current_directory = os.getcwd()
    os.chdir(WORK_PATH)
    try:
        subprocess.run(cmd, check=True)
    finally:
        os.chdir(current_directory)

def convert_to_pdb(input_file, output_file):
    _, ext = os.path.splitext(input_file)
    
    ext_lower = ext.lower()
    if ext_lower in ['.smi', '.mol2', '.cif', '.sdf']:
        input_format = ext_lower[1:]
    else:
        raise ValueError(f"Unsupported input format for conversion: {ext}")
    cmd = ["obabel", "-i", input_format, input_file, "-o", "pdb", "-O", output_file, "--gen3d"]
    subprocess.run(cmd, check=True)


def read_json(json_file, data_path=DATA_PATH, work_path=WORK_PATH):
    data = {}
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    receptor = os.path.join(data_path, data["receptor"])
    ligand = os.path.join(data_path, data["ligand"])
    output = os.path.join(work_path, data["output"])
    
    center_x = data["center"]["x"]
    center_y = data["center"]["y"]
    center_z = data["center"]["z"]
    
    size_x = data["size"]["x"]
    size_y = data["size"]["y"]
    size_z = data["size"]["z"]
    
    return receptor, ligand, output, center_x, center_y, center_z, size_x, size_y, size_z

def check_file_exists(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found!")
    return filename

def create_zip(files, zip_filename):
    with zipfile.ZipFile(zip_filename, 'w') as zipf:
        for file in files:
            zipf.write(file, os.path.basename(file))
    print(f"Saving to: {zip_filename}")


def run_docking(json_file, data_path=DATA_PATH, work_path=WORK_PATH, output_folder="output"):
    receptor, ligand, output, center_x, center_y, center_z, size_x, size_y, size_z = read_json(json_file)

    receptor = check_file_exists(receptor)
    ligand = check_file_exists(ligand)
    
    # Convert receptor if it's not in pdbqt format
    base, ext = os.path.splitext(receptor)
    if ext.lower() == '.pdb':
        # Convert receptor file to file with aminoacid sequence format only
        _, tail = os.path.split(receptor)
        receptor_aux = os.path.join(work_path, tail)
        shutil.copy2(receptor, receptor_aux)
        filtered_receptor = receptor_aux.replace(".pdb", "_filtered.pdb")
        filter_pdb(receptor_aux, filtered_receptor) 
        receptor = filtered_receptor

        receptor_pdbqt = os.path.abspath(base + '.pdbqt')
        prepare_receptor(receptor, receptor_pdbqt)
        check_file_exists(receptor_pdbqt)
        receptor = receptor_pdbqt

    # Convert ligand if it's not in pdbqt format
    ##### FIXED duplicities
    base, ext = os.path.splitext(ligand)
    _, tail = os.path.split(ligand)
    ext_lower = ext.lower()
    if ext_lower in ['.smi', '.mol2', '.cif', '.sdf', '.pdb']:
        ligand_pdb = os.path.join(WORK_PATH, tail.replace(ext, '.pdb'))
        if ext_lower != '.pdb':
            print(f"Converting {ligand} to {ligand_pdb} ...")
            convert_to_pdb(ligand, ligand_pdb)
            check_file_exists(ligand_pdb)
        ligand_pdbqt = os.path.join(WORK_PATH, tail.replace(ext, '.pdbqt'))
        #ligand_pdbqt = os.path.abspath(f'{WORK_DIR}{base}.pdbqt')
        print(f"Converting {ligand} to {ligand_pdbqt} ...")
        prepare_ligand(ligand_pdb, ligand_pdbqt)
        check_file_exists(ligand_pdbqt)
        ligand = ligand_pdbqt
    else:
        raise ValueError("Unsupported ligand file format: " + ext)
    
    # Create output folder
    os.makedirs(os.path.dirname(output), exist_ok=True)

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
    subprocess.run(vina_cmd, check=True)

    output_folder_full = os.path.join(data_path, output_folder)
    os.makedirs(output_folder_full, exist_ok=True)
    zip_filename = os.path.join(output_folder_full, "results.zip")
    files_to_zip = [output, receptor, ligand]  # ligand added to gzip
    create_zip(files_to_zip, zip_filename)

    print(f"Results saved to: {zip_filename}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 docking_script.py path_to_json_file")
        sys.exit(1)
    json_file = sys.argv[1]
    run_docking(json_file)
