from Bio.PDB import PDBParser, calc_dihedral
from Bio.PDB.Polypeptide import protein_letters_3to1, is_aa
import numpy as np
import os

def calculate_phi_psi_and_save_to_csv(pdb_file_path, output_csv_path):
    parser = PDBParser(QUIET=True)
    structure_id = os.path.basename(pdb_file_path).split('.')[0]
    structure = parser.get_structure(structure_id, pdb_file_path)
    
    # Open the output CSV file. 'w' mode will overwrite existing file, writing the header.
    with open(output_csv_path, 'w') as csv_file:
        # Write header
        csv_file.write("ResidueNumber,ResidueName,Phi,Psi\n")
        
        for model in structure:
            for chain in model:
                residues = list(chain)
                for i, res in enumerate(residues):
                    if not is_aa(res, standard=True):
                        continue

                    try:
                        phi = psi = None
                        if i > 0:
                            phi = calc_dihedral(residues[i-1]['C'].get_vector(),
                                                res['N'].get_vector(),
                                                res['CA'].get_vector(),
                                                res['C'].get_vector())

                        if i < len(residues) - 1:
                            psi = calc_dihedral(res['N'].get_vector(),
                                                res['CA'].get_vector(),
                                                res['C'].get_vector(),
                                                residues[i+1]['N'].get_vector())

                        if phi is not None and psi is not None:
                            res_num = res.get_id()[1]
                            try:
                                # Using protein_letters_3to1 dictionary
                                res_name = protein_letters_3to1[res.get_resname()]
                                # Write data directly to CSV
                                csv_file.write(f"{res_num},{res_name},{np.degrees(phi)},{np.degrees(psi)}\n")
                            except KeyError:
                                continue

                    except KeyError:
                        continue

def process_structures_and_save_angles(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for pdb_file in os.listdir(input_dir):
        if pdb_file.endswith(".pdb"):
            pdb_file_path = os.path.join(input_dir, pdb_file)
            print(f"Processing {pdb_file}...")
            output_file = pdb_file.replace('.pdb', '_angles.csv')
            output_csv_path = os.path.join(output_dir, output_file)
            calculate_phi_psi_and_save_to_csv(pdb_file_path, output_csv_path)
            print(f"Saved {output_file}")

# Ensure to set the correct paths for your input_directory and output_directory
input_directory = ""
output_directory = ""


process_structures_and_save_angles(input_directory, output_directory)
