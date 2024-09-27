
import os
import pandas as pd

def extract_asa_from_dssp(dssp_file, chain_id='A'):
    asa_data = []
    with open(dssp_file, 'r') as f:
        for line in f:
            columns = line.split()
            if len(columns) >= 5 and columns[0] == chain_id:
                resnum = int(columns[1].strip())
                res_type = columns[2].strip()
                asa_value = columns[4].strip()
                
                # Handle 'NA' or other invalid ASA values
                try:
                    asa = float(asa_value)
                except ValueError:
                    asa = 'NA'  # Or use 0.0, None, etc.
                    print(f"Residue {resnum} has an ASA value of 'NA' in file {dssp_file}")

                asa_data.append([resnum, res_type, asa])
    
    return pd.DataFrame(asa_data, columns=['Residue_Number', 'Residue_Type', 'ASA'])

def process_dssp_files(dssp_dir, output_dir, chain_id='A'):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for file in os.listdir(dssp_dir):
        if file.endswith('.dssp'):
            print(f"Processing file: {file}")
            dssp_file_path = os.path.join(dssp_dir, file)
            asa_df = extract_asa_from_dssp(dssp_file_path, chain_id)
            output_file_path = os.path.join(output_dir, file.replace('.dssp', '_asa.csv'))
            asa_df.to_csv(output_file_path, index=False)
            print(f"ASA data saved to: {output_file_path}")

if __name__ == '__main__':
    dssp_dir = ""  # Path to the directory containing DSSP files
    output_dir = ""  # Path to the output directory for CSV files
    process_dssp_files(dssp_dir, output_dir, chain_id='A')
