import os
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1  # Correct import for three-letter to one-letter conversion
import warnings

def get_residue_numbers(pdb_file_path, chain_id):
    parser = PDBParser(QUIET=True)

    # Check if the PDB file exists before attempting to parse it
    if not os.path.exists(pdb_file_path) or pdb_file_path.startswith('.') or not pdb_file_path.endswith('.pdb'):
        print(f"Skipping: File not found or invalid: {pdb_file_path}")
        return []

    structure = parser.get_structure('', pdb_file_path)
    
    # Check if the structure contains models
    if len(structure) == 0:
        print(f"No models found in {pdb_file_path}")
        return []

    model = structure[0]  # Assuming the first model in the PDB file

    residue_numbers = []
    if chain_id in model:
        for residue in model[chain_id]:
            if is_aa(residue, standard=True):
                residue_numbers.append((residue.id[1], seq1(residue.resname)))
    else:
        print(f"Chain ID {chain_id} not found in {pdb_file_path}")
        return []

    return residue_numbers


def parse_alignment(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    alignment = {"seq1": "", "seq2": "", "markers": ""}
    capturing = False
    for line in lines:
        if "denotes residue pairs" in line or "denotes other aligned residues" in line:
            capturing = True
            continue
        if capturing:
            if not alignment["seq1"]:
                alignment["seq1"] = line.strip()
            elif not alignment["markers"]:
                alignment["markers"] = line.strip()
            elif not alignment["seq2"]:
                alignment["seq2"] = line.strip()
                break
    return alignment["seq1"], alignment["markers"], alignment["seq2"]

def extract_ids_from_filename(filename):
    # Check if '_vs_' is in the filename to ensure it's a comparison file
    if '_vs_' not in filename:
        raise ValueError(f"Filename '{filename}' does not contain expected '_vs_' pattern.")

    parts = filename.split('_vs_')
    
    try:
        # Attempt to split the first part to extract PDB ID and Chain ID
        pdb_id_1, chain_part_1 = parts[0].rsplit('_', 1)
        chain_id_1 = chain_part_1.replace('Chain', '')
    except ValueError:
        raise ValueError(f"Error extracting PDB ID and Chain ID from the first part: '{parts[0]}'.")

    try:
        # Attempt to split the second part to extract PDB ID and Chain ID
        # Remove file extension before splitting
        pdb_id_2_chain_part_2 = parts[1].rsplit('.', 1)[0]
        pdb_id_2, chain_part_2 = pdb_id_2_chain_part_2.rsplit('_', 1)
        chain_id_2 = chain_part_2.replace('Chain', '')
    except ValueError:
        raise ValueError(f"Error extracting PDB ID and Chain ID from the second part: '{parts[1]}'.")

    return pdb_id_1, chain_id_1, pdb_id_2, chain_id_2


def compare_and_save_mappings(file_path, pdb_dir, output_dir):
    seq1, markers, seq2 = parse_alignment(file_path)
    pdb_id_1, chain_id_1, pdb_id_2, chain_id_2 = extract_ids_from_filename(os.path.basename(file_path))

    residue_numbers_1 = get_residue_numbers(os.path.join(pdb_dir, f"{pdb_id_1}_Chain{chain_id_1}.pdb"), chain_id_1)
    residue_numbers_2 = get_residue_numbers(os.path.join(pdb_dir, f"{pdb_id_2}_Chain{chain_id_2}.pdb"), chain_id_2)

    comparison_data = []
    index_1, index_2 = 0, 0

    for alignment_char, aa1, aa2 in zip(markers, seq1, seq2):
        if aa1 != '-':
            if index_1 < len(residue_numbers_1):
                res_num_1, aa_code_1 = residue_numbers_1[index_1]
                index_1 += 1
            else:
                print(f"Index out of range for sequence 1 at position {index_1}.")
                res_num_1, aa_code_1 = '-', '-'
        else:
            res_num_1, aa_code_1 = '-', '-'

        if aa2 != '-':
            if index_2 < len(residue_numbers_2):
                res_num_2, aa_code_2 = residue_numbers_2[index_2]
                index_2 += 1
            else:
                print(f"Index out of range for sequence 2 at position {index_2}.")
                res_num_2, aa_code_2 = '-', '-'
        else:
            res_num_2, aa_code_2 = '-', '-'

        if alignment_char != ' ':
            comparison_data.append(f"{res_num_1}{aa_code_1} - {res_num_2}{aa_code_2}")

    output_filename = os.path.splitext(os.path.basename(file_path))[0] + "_residue_mapping.txt"
    save_residue_mapping(output_dir, output_filename, comparison_data)


def save_residue_mapping(output_dir, filename, comparison_data):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, filename)
    with open(output_path, 'w') as f:
        for line in comparison_data:
            f.write(line + "\n")


# Example usage
pdb_dir = ""
tm_align_output_dir = ""
output_dir = ""

for filename in os.listdir(tm_align_output_dir):
    # Skip hidden or non-TXT files
    if filename.startswith('.') or not filename.endswith(".txt"):
        print(f"Skipping hidden or non-TXT file: {filename}")
        continue

    file_path = os.path.join(tm_align_output_dir, filename)
    compare_and_save_mappings(file_path, pdb_dir, output_dir)
    print(f"Residue mapping for {filename} saved.")
