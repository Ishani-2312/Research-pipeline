import os
import pandas as pd
import warnings
from Bio.PDB import PDBParser, DSSP

# Suppress specific warnings
warnings.filterwarnings('ignore', message="Ignoring unrecognized record 'TER'")
warnings.filterwarnings('ignore', message="DSSP could not be created due to an error")

def run_dssp(pdb_file, dssp_executable):
    p = PDBParser(QUIET=True)
    try:
        structure = p.get_structure('X', pdb_file)
        model = structure[0]  # Assuming first model
        dssp = DSSP(model, pdb_file, dssp=dssp_executable)
        return dssp
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return None

def parse_dssp_data(dssp):
    ss_map = {
        'H': 'Helix', 'B': 'Beta Strand', 'E': 'Beta Strand', 
        'G': 'Helix', 'I': 'Helix', 'T': 'Turn', 'S': 'Coil', ' ': 'Coil'
    }
    ss_counts = {'Helix': 0, 'Beta Strand': 0, 'Turn': 0, 'Coil': 0}
    results = []
    
    for key in dssp.keys():
        aa = dssp[key][1]
        ss = ss_map.get(dssp[key][2], 'Coil')
        resnum = key[1][1]
        results.append({'ResidueNumber': resnum, 'AminoAcid': aa, 'SecondaryStructure': ss})
        ss_counts[ss] += 1

    total_residues = sum(ss_counts.values())
    ss_percentages = {ss: (count / total_residues) * 100 for ss, count in ss_counts.items()}
    
    return results, ss_percentages

def process_pdb_files(input_directory, output_directory, dssp_executable):
    os.makedirs(output_directory, exist_ok=True)

    for filename in os.listdir(input_directory):
        if filename.endswith(".pdb"):
            pdb_file_path = os.path.join(input_directory, filename)
            protein_name = os.path.splitext(filename)[0]
            output_file = os.path.join(output_directory, f"SS_{protein_name}.csv")
            
            dssp = run_dssp(pdb_file_path, dssp_executable)
            if dssp is not None:
                results, ss_percentages = parse_dssp_data(dssp)
                
                # Save results to CSV
                df_results = pd.DataFrame(results)
                df_ss_percentages = pd.DataFrame.from_dict(ss_percentages, orient='index', columns=['Percentage']).reset_index()
                df_ss_percentages = df_ss_percentages.rename(columns={'index': 'SecondaryStructure'})
                
                # Writing both results and percentages into the CSV
                with open(output_file, 'w') as f:
                    df_results.to_csv(f, index=False)
                    f.write("\n\nSecondary Structure Percentages\n")
                    df_ss_percentages.to_csv(f, index=False)
                
                print(f"DSSP analysis completed for {protein_name}. Results saved to {output_file}.")
            else:
                print(f"Failed to process {pdb_file_path}: DSSP failed to produce an output.")

if __name__ == "__main__":
    input_directory = "/mnt/"  # Update with your PDB files directory
    output_directory = "/mnt/"    # Update with desired output directory
    dssp_executable = "/usr/local/bin/mkdssp"  # Update with the full path to your DSSP executable
    process_pdb_files(input_directory, output_directory, dssp_executable)
