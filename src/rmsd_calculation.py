import os
import pandas as pd
import numpy as np
from scipy.spatial import distance

# Directory containing the Excel coordinate files
coord_dir = ""  # Replace with your directory path
output_dir = ""  # Replace with your desired output directory

# Ensure the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Function to sort residues
def sort_residues(residues):
    def alphanum_key(key):
        import re
        return [int(text) if text.isdigit() else text for text in re.split('([0-9]+)', key)]
    return sorted(residues, key=alphanum_key)

# Function to process each Excel file
def process_excel_file(file_path):
    # Read the Excel file
    df = pd.read_excel(file_path)
    
    # Identify the reference model (assuming it's the first model in the file)
    reference_model = df['model'].iloc[0]
    print(f"Using {reference_model} as the reference model")

    # Define the atoms to consider for distance calculations
    atoms_to_consider = ['CA']

    # Create a DataFrame to store the results
    result_df = pd.DataFrame(columns=['Residue', 'Amino_Acid'] + [f"RMSD:{model}" for model in df['model'].unique() if model != reference_model])

    # Filter rows for the reference structure
    reference_df = df[df['model'] == reference_model]

    # Extract all unique models excluding the reference
    unique_models = df['model'].unique()
    unique_models = [model for model in unique_models if model != reference_model]

    # Iterate over unique models
    for model in unique_models:
        current_df = df[df['model'] == model]
        for residue in reference_df['resi'].unique():
            ref_atoms_df = reference_df[(reference_df['resi'] == residue) & (reference_df['name'].isin(atoms_to_consider))]
            curr_atoms_df = current_df[(current_df['resi'] == residue) & (current_df['name'].isin(atoms_to_consider))]

            if not ref_atoms_df.empty and not curr_atoms_df.empty:
                ref_coords = ref_atoms_df[['x', 'y', 'z']].values
                curr_coords = curr_atoms_df[['x', 'y', 'z']].values
                distances = distance.cdist(ref_coords, curr_coords, 'euclidean')

                # Compute RMSD
                rmsd = np.sqrt(np.mean(np.square(distances)))
                
                # Update or insert the row
                if str(residue) in result_df['Residue'].values:
                    result_df.loc[result_df['Residue'] == str(residue), f"RMSD:{model}"] = rmsd
                else:
                    new_row = pd.DataFrame({
                        'Residue': [str(residue)],
                        'Amino_Acid': [ref_atoms_df['resn'].iloc[0]],
                        f"RMSD:{model}": [rmsd]
                    })
                    if not new_row.isna().all().all():
                        result_df = pd.concat([result_df, new_row], ignore_index=True)

    # Sort the DataFrame by Residue using the custom sorting function
    sorted_residues = sort_residues(result_df['Residue'].unique())
    result_df['Residue'] = pd.Categorical(result_df['Residue'], categories=sorted_residues, ordered=True)
    result_df = result_df.sort_values(by='Residue').reset_index(drop=True)

    # Calculate the mean RMSD for each residue across all models
    result_df['Mean_RMSD'] = result_df.loc[:, [col for col in result_df.columns if col.startswith('RMSD:')]].mean(axis=1)

    # Save the results to a CSV file
    protein_name = os.path.splitext(os.path.basename(file_path))[0].replace('_coord', '')
    output_file = os.path.join(output_dir, f"results_mean_RMSD_{protein_name}.csv")
    result_df.to_csv(output_file, index=False)

    print(f"RMSD calculation completed for {protein_name} and results saved to {output_file}")

# Process all Excel files in the directory
for file in os.listdir(coord_dir):
    if file.endswith("_coord.xlsx"):
        process_excel_file(os.path.join(coord_dir, file))
