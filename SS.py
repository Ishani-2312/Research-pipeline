import os
import pandas as pd
from Bio.PDB import PDBParser, DSSP

def run_dssp(pdb_file, dssp_executable):
    p = PDBParser()
    structure = p.get_structure('X', pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp=dssp_executable)
    return dssp

def parse_dssp_data(dssp):
    ss_map = {'H': 'Helix', 'B': 'Beta Strand', 'E': 'Beta Strand', 'G': 'Helix', 'I': 'Helix', 'T': 'Turn', 'S': 'Coil', ' ': 'Coil'}
    ss_counts = {'Helix': 0, 'Beta Strand': 0, 'Turn': 0, 'Coil': 0}
    results = []
    
    for key in dssp.keys():
        aa = dssp[key][1]
        ss = ss_map.get(dssp[key][2], 'Coil')
        resnum = key[1][1]
        results.append({'ResidueNumber': resnum, 'AminoAcid': aa, 'SecondaryStructure': ss})
        ss_counts[ss] += 1

    return results, ss_counts

def calculate_ss_percentages(ss_counts, total_residues):
    return {ss: (count / total_residues) * 100 for ss, count in ss_counts.items()}

def main(pdb_file, output_file, dssp_executable):
    dssp = run_dssp(pdb_file, dssp_executable)
    results, ss_counts = parse_dssp_data(dssp)
    
    total_residues = sum(ss_counts.values())
    ss_percentages = calculate_ss_percentages(ss_counts, total_residues)
    
    df_results = pd.DataFrame(results)
    
    # Create a DataFrame for the percentages row
    percentages_row = pd.DataFrame([{
        'ResidueNumber': 'Total', 
        'AminoAcid': '', 
        'SecondaryStructure': f"Helix: {ss_percentages['Helix']:.2f}%, Beta Strand: {ss_percentages['Beta Strand']:.2f}%, Turn: {ss_percentages['Turn']:.2f}%, Coil: {ss_percentages['Coil']:.2f}%"
    }])
    
    # Concatenate the percentages row to the results DataFrame
    df_results = pd.concat([df_results, percentages_row], ignore_index=True)
    
    df_results.to_csv(output_file, index=False)

    print(f"DSSP analysis completed. Results saved to {output_file}.")
    print(f"Secondary Structure Percentages: {ss_percentages}")

if __name__ == "__main__":
    pdb_file = "/mnt/c/Users/msc2/Desktop/Project_data/mhc1/pdb_chainA/6O53_ChainA.pdb"  # Update with your PDB file path
    output_file = "/mnt/c/Users/msc2/Desktop/Project_data/mhc1/ss_output.csv"  # Update with desired output CSV file path
    dssp_executable = "/usr/local/bin/mkdssp"  # Update with the full path to your DSSP executable
    main(pdb_file, output_file, dssp_executable)
