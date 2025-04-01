import os
import subprocess
import json
import sys
import zipfile
import shutil
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select
from rdkit import Chem
from rdkit.Chem import AllChem

MGLTOOLS_LIB_PATH = '/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs'
MGLTOOLS_BIN_PATH = '/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/'
DATA_PATH = "/data/"
WORK_PATH = "/tmp/"

PREPARATION_TIMEOUT_SECONDS = 600
DOCKING_TIMEOUT_SECONDS = 900

def is_valid_pdbqt(pdbqt_file):
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    if abs(x) > 1e-3:
                        return True
                except Exception:
                    continue
    return False

class MySelect(Select):
    def accept_residue(self, residue):
        return residue.get_resname() in [
            "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS",
            "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",
            "TRP", "TYR", "HID", "HSP", "HIE", "HIP", "CYX", "CSS"
        ]

def filter_cif_file_biopython(input_cif, output_cif):
    # Filter the CIF file to keep only selected sections
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
        # Skip the first two lines (data_ID and #)
        for _ in range(2):
            next(infile)

        outfile.write("data_ID\n# \n")

        for key, value in dict_saved_lines.items():
            for line in value:
                outfile.write(line + "\n")
            outfile.write("# \n")
        
        outfile.writelines(infile)
    
    os.remove(temp_cif)

def filter_pdb(input_pdb, output_pdb):
    # Filter the PDB file to keep only ATOM records
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        writing_hetatm = False
        for line in f_in:
            if line.startswith("ATOM"):
                f_out.write(line)
                writing_hetatm = True 
            elif line.startswith("HETATM") and writing_hetatm:
                modified_line = "ATOM  " + line[6:]
                f_out.write(modified_line)
            elif line.startswith("TER") or line.startswith("END"):
                writing_hetatm = False
                break 
    print(f"Filtered PDB for inspection saved to {output_pdb}")

def prepare_receptor(pdb_file, pdbqt_file, use_timeout):
    # Prepare the receptor using MGLTools prepare_receptor4.py
    cmd = ["/usr/local/bin/python2.7",
           os.path.join(MGLTOOLS_BIN_PATH, 'prepare_receptor4.py'),
           '-A', "hydrogens", '-U', "nphs _lps_waters_deleteAltB",
           '-r', pdb_file, '-o', pdbqt_file]

    if use_timeout:
        subprocess.run(cmd, check=True, timeout=PREPARATION_TIMEOUT_SECONDS)
    else:
        subprocess.run(cmd, check=True)

def prepare_ligand(input_file, output_file, use_timeout):
    # Prepare the ligand using MGLTools prepare_ligand4.py
    cmd = ["/usr/local/bin/python2.7", os.path.join(MGLTOOLS_BIN_PATH, 'prepare_ligand4.py'), '-l', input_file, '-o', output_file]
    current_directory = os.getcwd()
    os.chdir(WORK_PATH)
    try:
        if use_timeout:
            subprocess.run(cmd, check=True, timeout=PREPARATION_TIMEOUT_SECONDS)
        else:
            subprocess.run(cmd, check=True)
    finally:
        os.chdir(current_directory)

def convert_to_pdb(input_file, output_file, use_timeout):
    # Convert ligand to PDB using Open Babel or RDKit if input is SMILES
    _, ext = os.path.splitext(input_file)
    
    ext_lower = ext.lower()
    if ext_lower == ".smi":
        with open(input_file, 'r') as f:
            smiles = f.read().strip()

        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            raise ValueError("Invalid SMILES string.")
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception as e:
            print("Warning: Kekulization failed:", e)

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)
        pdb_block = Chem.MolToPDBBlock(mol)

        with open(output_file, 'w') as f:
            f.write(pdb_block)

    elif ext_lower in ['.mol2', '.cif', '.sdf', '.pdb']:
        input_format = ext_lower[1:]
        cmd = ["obabel", "-i", input_format, input_file,
               "-o", "pdb", "-O", output_file, "--gen3d"]
        if use_timeout:
            subprocess.run(cmd, check=True, timeout=PREPARATION_TIMEOUT_SECONDS)
        else:
            subprocess.run(cmd, check=True)
    else:
        raise ValueError(f"Unsupported input format: {ext}")

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

# ---------------------------
# Fallback mechanism using RDKit (with multiconformer embedding)
# ---------------------------
def add_gasteiger_charges_with_rdkit(mol):
    Chem.SanitizeMol(mol)
    AllChem.ComputeGasteigerCharges(mol)
    return mol

def fallback_rdkit_conversion(original_smiles_file):
    # Read the SMILES string from the file and convert using RDKit
    with open(original_smiles_file, 'r', encoding='utf-8-sig') as f:
        content = f.read().strip()

    if not content:
        raise ValueError("Fallback: File is empty.")

    candidate = content.splitlines()[0].strip()
    print(f"Extracted candidate SMILES: '{candidate}'")

    mol = Chem.MolFromSmiles(candidate)
    if mol is None:
        raise ValueError(f"Fallback: Invalid SMILES string: '{candidate}'")
    mol = Chem.AddHs(mol)

    # Embed multiple conformers using ETKDG
    params = AllChem.ETKDG()
    params.randomSeed = 42

    cids = AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params)
    if not cids:
        raise Exception("Fallback: EmbedMultipleConfs failed.")

    energies = AllChem.UFFOptimizeMoleculeConfs(mol)
    best_cid = min(range(len(energies)), key=lambda i: energies[i])

    conf = mol.GetConformer(best_cid)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)

    AllChem.UFFOptimizeMolecule(mol)
    AllChem.ComputeGasteigerCharges(mol)

    return mol

def rdkit_mol_to_pdbqt(mol, output_file, resname="UNL", chain="A", resid=1):
    # Convert RDKit molecule to PDBQT format compatible with AutoDock Vina
    pdbqt_lines = []
    pdbqt_lines.append("REMARK   Generated by rdkit_mol_to_pdbqt")
    pdbqt_lines.append("ROOT")

    for i, atom in enumerate(mol.GetAtoms(), start=1):
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        try:
            charge = float(atom.GetProp("_GasteigerCharge"))
        except Exception:
            charge = 0.0
        element = atom.GetSymbol()
        atom_name = element
        line = f"ATOM  {i:5d}  {atom_name:<4s}{resname:>3s} {chain}{resid:4d}    {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00 {charge:6.3f}          {element:>2s}"
        pdbqt_lines.append(line)

    pdbqt_lines.append("ENDROOT")
    pdbqt_lines.append("TORSDOF 0")
    pdbqt_lines.append("END")

    with open(output_file, 'w') as f:
        f.write("\n".join(pdbqt_lines))

    print(f"Custom RDKit -> PDBQT conversion complete, file saved to {output_file}")

def prepare_ligand_with_fallback(input_file, original_smiles_file, output_file):
    try:
        prepare_ligand(input_file, output_file, use_timeout=False)
        if not is_valid_pdbqt(output_file):
            raise RuntimeError("PDBQT file generated by MGLTools contains invalid coordinates.")

    except Exception as e:
        print(f"MGLTools ligand preparation failed: {e}")
        print("Running fallback using RDKit...")
        mol = fallback_rdkit_conversion(original_smiles_file)

        rdkit_mol_to_pdbqt(mol, output_file)
        if not is_valid_pdbqt(output_file):
            raise RuntimeError("Fallback: Custom PDBQT conversion is invalid.")

        print("Fallback conversion was successful.")

# ---------------------------
# Main workflow
# ---------------------------
def run_docking(json_file, data_path=DATA_PATH, work_path=WORK_PATH, output_folder="output", use_timeout=True):
    receptor, ligand, output, center_x, center_y, center_z, size_x, size_y, size_z = read_json(json_file)
    receptor = check_file_exists(receptor)
    ligand = check_file_exists(ligand)
    
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
    prepare_receptor(receptor, receptor_pdbqt, use_timeout)
    check_file_exists(receptor_pdbqt)
    receptor = receptor_pdbqt

    base, ext = os.path.splitext(ligand)
    _, tail = os.path.split(ligand)
    ext_lower = ext.lower()

    if ext_lower in ['.smi', '.mol2', '.cif', '.sdf', '.pdb']:
        ligand_pdb = os.path.join(work_path, tail.replace(ext, '.pdb'))
        if ext_lower != '.pdb':
            print(f"Converting {ligand} to {ligand_pdb} ...")
            convert_to_pdb(ligand, ligand_pdb, use_timeout)
            check_file_exists(ligand_pdb)

        ligand_pdbqt = os.path.join(work_path, tail.replace(ext, '.pdbqt'))
        print(f"Converting {ligand} to {ligand_pdbqt} ...")

        try:
            prepare_ligand(ligand_pdb, ligand_pdbqt, use_timeout)
            if not is_valid_pdbqt(ligand_pdbqt):
                raise RuntimeError("Invalid coordinates in PDBQT generated by MGLTools.")
            check_file_exists(ligand_pdbqt)

        except Exception as e:
            print(f"MGLTools ligand preparation failed: {e}")
            print("Running fallback using RDKit...")
            mol = fallback_rdkit_conversion(ligand)
            rdkit_mol_to_pdbqt(mol, ligand_pdbqt)
            if not is_valid_pdbqt(ligand_pdbqt):
                raise RuntimeError("Fallback: Custom PDBQT conversion is invalid.")
            check_file_exists(ligand_pdbqt)

        ligand = ligand_pdbqt
    else:
        raise ValueError("Unsupported ligand file format: " + ext)

    os.makedirs(os.path.dirname(output), exist_ok=True)

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
    subprocess.run(vina_cmd, check=True, timeout=DOCKING_TIMEOUT_SECONDS)

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
    run_docking(json_file, use_timeout=False)
