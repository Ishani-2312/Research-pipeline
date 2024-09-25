import os
import pandas as pd
from Bio.PDB import PDBParser, PDBExceptions

def extract_b_factors(pdb_file):
    """Extract average B-factors for each residue from a PDB file."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('structure', pdb_file)
    except (PDBExceptions.PDBConstructionException, FileNotFoundError) as e:
        print(f"Error reading PDB file {pdb_file}: {e}")
        return []
    
    b_factors = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                res_b_factors = [atom.get_bfactor() for atom in residue if atom.is_disordered() == 0]
                if res_b_factors:
                    avg_b_factor = sum(res_b_factors) / len(res_b_factors)
                    b_factors.append({
                        'Model': model.id,
                        'Chain': chain.id,
                        'ResidueName': residue.resname,
                        'ResidueNumber': residue.id[1],
                        'AverageBFactor': avg_b_factor
                    })
    
    return b_factors

def process_pdb_files(input_dir, output_dir):
    """Process all PDB files in the input directory and save B-factors to the output directory."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for file in os.listdir(input_dir):
        if file.endswith('.pdb'):
            pdb_file = os.path.join(input_dir, file)
            b_factors = extract_b_factors(pdb_file)
            
            if b_factors:
                df_b_factors = pd.DataFrame(b_factors)
                output_file = os.path.join(output_dir, f"{os.path.splitext(file)[0]}.csv")
                df_b_factors.to_csv(output_file, index=False)
                print(f"B-factors extracted and saved to {output_file}")
            else:
                print(f"No B-factors found or error in file: {file}")

if __name__ == "__main__":
    input_dir = "C:\\Users\\msc2\\Desktop\\Project_data\\ENTROPY\\pdb_ref"  # Update with your input directory path
    output_dir = "C:\\Users\\msc2\\Desktop\\Project_data\\ENTROPY\\b-factor"  # Update with your desired output directory path
    process_pdb_files(input_dir, output_dir)
