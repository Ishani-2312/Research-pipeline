import os
import pandas as pd
import numpy as np
from scipy.stats import entropy
import matplotlib.pyplot as plt
import seaborn as sns
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_histogram_entropy(data, bins):
    """
    Calculate the entropy of a histogram of the given data.

    Parameters:
    - data (Series): The data for which to calculate the entropy.
    - bins (int): The number of bins for the histogram.

    Returns:
    - float: The calculated entropy.
    """
    data = pd.to_numeric(data, errors='coerce').dropna()  # Ensure data is numeric and drop NaNs
    if data.empty:
        logging.warning("Empty data series provided for entropy calculation.")
        return np.nan  # Return NaN for empty data

    histogram, bin_edges = np.histogram(data, bins=bins, range=(0, 360), density=True)
    histogram = histogram[histogram > 0]

    return entropy(histogram, base=2)  # Base 2 logarithm for binary entropy

def plot_kde(data_phi, data_psi, residue_number, output_dir, protein_name):
    """
    Plot the KDE for both Δϕ and Δψ for a given residue.

    Parameters:
    - data_phi (Series): Δϕ data for the residue.
    - data_psi (Series): Δψ data for the residue.
    - residue_number (int): The residue number for labeling plots.
    - output_dir (str): The directory to save the plots.
    - protein_name (str): The name of the protein for labeling plots.
    """
    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    sns.histplot(data_phi, kde=True, color='blue')
    plt.title(f'KDE of Δϕ for Residue {residue_number}')
    plt.xlabel('Δϕ (Degrees)')
    plt.ylabel('Density')
    
    plt.subplot(1, 2, 2)
    sns.histplot(data_psi, kde=True, color='blue')
    plt.title(f'KDE of Δψ for Residue {residue_number}')
    plt.xlabel('Δψ (Degrees)')
    plt.ylabel('Density')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{protein_name}_Residue_{residue_number}_KDE.png"))
    plt.close()

def plot_histogram(data_phi, data_psi, residue_number, output_dir, protein_name, bins):
    """
    Plot the histogram for both Δϕ and Δψ for a given residue.

    Parameters:
    - data_phi (Series): Δϕ data for the residue.
    - data_psi (Series): Δψ data for the residue.
    - residue_number (int): The residue number for labeling plots.
    - output_dir (str): The directory to save the plots.
    - protein_name (str): The name of the protein for labeling plots.
    - bins (int): The number of bins for the histogram.
    """
    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    plt.hist(data_phi, bins=bins, range=(0, 360), color='purple', alpha=0.7, edgecolor='black')
    plt.title(f'Histogram of Δϕ for Residue {residue_number}')
    plt.xlabel('Δϕ (Degrees)')
    plt.ylabel('Density')
    
    plt.subplot(1, 2, 2)
    plt.hist(data_psi, bins=bins, range=(0, 360), color='purple', alpha=0.7, edgecolor='black')
    plt.title(f'Histogram of Δψ for Residue {residue_number}')
    plt.xlabel('Δψ (Degrees)')
    plt.ylabel('Density')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{protein_name}_Residue_{residue_number}_Histogram.png"))
    plt.close()

def process_files(input_dir, output_dir, max_residue_number_dict, bins):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        if file_name.startswith("agg_") and file_name.endswith(".csv"):
            protein_name = file_name.split('_')[1].replace('.csv', '')
            input_file = os.path.join(input_dir, file_name)
            output_file = os.path.join(output_dir, f"entropy_{protein_name}.csv")

            max_residue_number = max_residue_number_dict.get(protein_name, None)

            # Load data
            try:
                data = pd.read_csv(input_file)
                logging.info(f"Processing file: {input_file}")

                # Filter residues up to max_residue_number if specified
                if max_residue_number is not None:
                    data = data[data['RefResidueNumber'] <= max_residue_number]

                # Ensure DeltaPhi and DeltaPsi columns are numeric
                data['DeltaPhi'] = pd.to_numeric(data['DeltaPhi'], errors='coerce')
                data['DeltaPsi'] = pd.to_numeric(data['DeltaPsi'], errors='coerce')

                # Calculate entropy for each residue number
                phi_entropy = data.groupby('RefResidueNumber')['DeltaPhi'].apply(lambda x: calculate_histogram_entropy(x, bins=bins)).reset_index()
                psi_entropy = data.groupby('RefResidueNumber')['DeltaPsi'].apply(lambda x: calculate_histogram_entropy(x, bins=bins)).reset_index()

                # Merge entropy results
                entropy_df = pd.merge(phi_entropy, psi_entropy, on='RefResidueNumber', suffixes=('_PhiEntropy', '_PsiEntropy'))
                entropy_df.columns = ['RefResidueNumber', 'PhiEntropy', 'PsiEntropy']

                # Save results, sorting by RefResidueNumber
                entropy_df.to_csv(output_file, index=False)
                logging.info(f"Entropy calculations completed and saved to {output_file}")

                # Plot distributions
                for residue_number in data['RefResidueNumber'].unique():
                    residue_data = data[data['RefResidueNumber'] == residue_number]
                    plot_kde(residue_data['DeltaPhi'], residue_data['DeltaPsi'], residue_number, output_dir, protein_name)
                    plot_histogram(residue_data['DeltaPhi'], residue_data['DeltaPsi'], residue_number, output_dir, protein_name, bins=bins)

            except Exception as e:
                logging.error(f"An error occurred while processing {input_file}: {e}")

if __name__ == "__main__":
    input_dir = "Z:\\msc2\\Desktop\\Project_data\\ENTROPY\\aggregate"
    output_dir = "Z:\\msc2\\Desktop\\Project_data\\ENTROPY\\Entropy3_fixed-60"
    
    # Dictionary specifying max residue numbers for specific protein families
    max_residue_number_dict = {
        'Carbonic-Anhydrase': 227,
        'Beta-Lactamase': 264,
        'Cytochrome-P450': 400,
        'Glutathione-S-Transferase': 200,
        'MHC1': 274,
        'Superoxide-Dismutase': 201,
        'Thioredoxin': 139,
        'Trypsin': 245,
        'Peroxidase':294,
        'Myoglobin':151,
        'p53':290,
        'Cytochrome-C':104,
        'KRAS':166,
        'Lysozyme-C':146,
        
        # Add other protein families and their max residue numbers here
    }

    process_files(input_dir, output_dir, max_residue_number_dict, bins=60)
