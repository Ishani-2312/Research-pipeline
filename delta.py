import pandas as pd
import os
import glob

def calculate_circular_difference(angle1, angle2):
    difference = (angle1 - angle2 + 180) % 360 - 180
    return abs(difference)

def load_phi_psi_values(file_path):
    df = pd.read_csv(file_path)
    df.set_index('ResidueNumber', inplace=True)
    return df

def parse_alignment_file(aln_file):
    alignment_mapping = {}
    with open(aln_file, 'r') as file:
        for line in file:
            parts = line.strip().split(' - ')
            # Check if the line contains residue mappings without gaps
            if '-' not in parts[0] and '-' not in parts[1]:
                # Extract the number, assuming the format 'numberLetter'
                ref_index = int(parts[0][:-1])
                target_index = int(parts[1][:-1])
                alignment_mapping[ref_index] = target_index
    return alignment_mapping


def calculate_deltas_and_save(aln_file, ref_df, target_df, output_file_path):
    delta_phi_psi = []
    alignment_mapping = parse_alignment_file(aln_file)
    for ref_index, target_index in alignment_mapping.items():
        if ref_index in ref_df.index and target_index in target_df.index:
            delta_phi = calculate_circular_difference(ref_df.at[ref_index, 'Phi'], target_df.at[target_index, 'Phi'])
            delta_psi = calculate_circular_difference(ref_df.at[ref_index, 'Psi'], target_df.at[target_index, 'Psi'])
            delta_phi_psi.append([ref_index, alignment_mapping[ref_index], delta_phi, delta_psi])
    delta_df = pd.DataFrame(delta_phi_psi, columns=['RefResidueNumber', 'TargetResidueNumber', 'DeltaPhi', 'DeltaPsi'])
    delta_df.to_csv(output_file_path, index=False)

def process_alignment_pairs(alignment_files_dir, phi_psi_files_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for aln_file in glob.glob(os.path.join(alignment_files_dir, '*_residue_mapping.txt')):
        base_name = os.path.basename(aln_file).replace('_residue_mapping.txt', '')
        ref_id, target_id = base_name.split('_vs_')
        
        ref_phi_psi_file = f"{ref_id}_angles.csv"
        target_phi_psi_file = f"{target_id}_angles.csv"
        
        ref_phi_psi_path = os.path.join(phi_psi_files_dir, ref_phi_psi_file)
        target_phi_psi_path = os.path.join(phi_psi_files_dir, target_phi_psi_file)
        
        if os.path.exists(ref_phi_psi_path) and os.path.exists(target_phi_psi_path):
            ref_df = load_phi_psi_values(ref_phi_psi_path)
            target_df = load_phi_psi_values(target_phi_psi_path)
            output_file_name = f"{ref_id}_vs_{target_id}_delta_phi_psi.csv"
            output_file_path = os.path.join(output_dir, output_file_name)
            calculate_deltas_and_save(aln_file, ref_df, target_df, output_file_path)
            print(f"Processed {output_file_name}")
        else:
            print(f"Missing phi/psi files for {base_name}")

# Set your directories here
alignment_files_dir = "C:\\Users\\msc2\\Desktop\\Project_data\\mhc1\\residue_mapping"
phi_psi_files_dir = "C:\\Users\\msc2\\Desktop\\Project_data\\mhc1\\dihedral"
output_dir = "C:\\Users\\msc2\\Desktop\\Project_data\\mhc1\\delta"

# Process each pair
process_alignment_pairs(alignment_files_dir, phi_psi_files_dir, output_dir)
