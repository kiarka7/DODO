import os
import subprocess
import json
import sys
import zipfile
import shutil
from Bio.PDB import MMCIFParser, Select, MMCIFIO

MGLTOOLS_LIB_PATH = '/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs'
MGLTOOLS_BIN_PATH = '/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/'
DATA_PATH = "/data/"
WORK_PATH ="/tmp/"

class MySelect(Select):
    def accept_residue(self, residue):
        return residue.get_resname() in ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

def filter_cif_file_biopython(input_cif, output_cif):
    # TODO: consider keeping STRUCT_CONN category record
    needed_sections = ['_entry.id', '_cell', '_symmetry']
    dict_saved_lines = {}

    with open(input_cif, 'r') as file:
        for line in file:
            for section in needed_sections:
                if line.startswith(section):
                    if dict_saved_lines.get(section) is None:
                        dict_saved_lines[section] = [line.strip()]
                    else:
                        dict_saved_lines[section].append(line.strip())

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('ID', input_cif)
    io = MMCIFIO()
    io.set_structure(structure)

    temp_cif = f"{output_cif}.tmp"

    io.save(temp_cif, select=MySelect())

    with open(temp_cif, 'r') as infile, open(output_cif, 'w') as outfile:
        # Skip the first two lines (to remove data_ID and #)
        for _ in range(2):
            next(infile)

        outfile.write("data_ID\n# \n")

        # Now, write the saved lines back to the top of the file
        for key, value in dict_saved_lines.items():
            for line in value:
                outfile.write(line + "\n")
            outfile.write("# \n")
        
        # Write the remaining lines to the new file
        outfile.writelines(infile)
    
    os.remove(temp_cif)

def filter_pdb(input_pdb, output_pdb):
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                f_out.write(line)

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
    
    # Convert receptor to file with aminoacid sequence format only
    base, ext = os.path.splitext(receptor)
    print(f"Cleaning receptor file {receptor} ...")

    if ext.lower() == ".cif":
        filtered_receptor = base + "_filtered.cif"
        filter_cif_file_biopython(receptor, filtered_receptor)
    elif ext.lower() == ".pdb":
        filtered_receptor = base + "_filtered.pdb"
        filter_pdb(receptor, filtered_receptor)
    else:
        raise ValueError("Unsupported receptor file format: " + ext)
    
    receptor = filtered_receptor
    print(f"Cleaned receptor saved to {receptor}")

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
