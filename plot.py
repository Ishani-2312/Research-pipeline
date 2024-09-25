import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def plot_entropy_and_ss(entropy_file, ss_file, output_dir, protein_name):
    try:
        # Load the entropy data
        entropy_data = pd.read_csv(entropy_file)

        # Load the secondary structure data, ignoring the last row if it contains summary information
        ss_data = pd.read_csv(ss_file)
        if isinstance(ss_data.iloc[-1, 0], str) and ss_data.iloc[-1, 0].startswith('Total'):
            ss_data = ss_data.iloc[:-1]

        # Ensure the columns to be merged have the same data type
        entropy_data['RefResidueNumber'] = entropy_data['RefResidueNumber'].astype(str)
        ss_data['ResidueNumber'] = ss_data['ResidueNumber'].astype(str)

        # Merge the two datasets on the residue number
        merged_data = pd.merge(entropy_data, ss_data, left_on='RefResidueNumber', right_on='ResidueNumber')

        # Convert merged residue numbers to numeric for plotting
        merged_data['RefResidueNumber'] = pd.to_numeric(merged_data['RefResidueNumber'], errors='coerce')

        # Define colors for secondary structures
        ss_colors = {'Helix': 'red', 'Beta Strand': 'blue', 'Coil': 'green', 'Turn': 'orange'}
        merged_data['SS_Color'] = merged_data['SecondaryStructure'].map(ss_colors)

        # Plotting
        fig, ax1 = plt.subplots(figsize=(16, 8))  # Adjusted to a more proportional size

        # Plot entropy lines with darker colors
        ax1.plot(merged_data['RefResidueNumber'], merged_data['PhiEntropy'], label='Phi Entropy', color='darkblue', linewidth=2)
        ax1.plot(merged_data['RefResidueNumber'], merged_data['PsiEntropy'], label='Psi Entropy', color='darkred', linewidth=2)

        # Plot secondary structure as color bars with lower opacity
        for i, row in merged_data.iterrows():
            ax1.axvspan(row['RefResidueNumber'] - 0.5, row['RefResidueNumber'] + 0.5, color=row['SS_Color'], alpha=0.2)

        # Set y-axis limits
        ax1.set_ylim(0, 5)  # Fixed y-axis from 0 to 8


        # Labels without title and legends
        ax1.set_xlabel('Residue Number', fontsize=24)  # Increased font size
        ax1.set_ylabel('Entropy', fontsize=24)  # Increased font size

        # Increasing font sizes for ticks
        ax1.tick_params(axis='both', which='major', labelsize=20)

        # Set x-ticks at regular intervals to avoid clutter
        interval = max(1, len(merged_data) // 20)  # Adjust the interval as needed
        xticks = merged_data['RefResidueNumber'][::interval]
        ax1.set_xticks(xticks)
        ax1.set_xticklabels(xticks, rotation=45, ha='right')

        # Add grid lines
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)

        plt.tight_layout()

        # Save the plot to a file with high quality
        output_file = os.path.join(output_dir, f"plot_{protein_name}.png")
        plt.savefig(output_file, dpi=600)

        # Show the plot (optional)
        plt.show()
        print(f"Plot generated and saved for {protein_name}")

    except Exception as e:
        print(f"An error occurred while processing {protein_name}: {e}")

if __name__ == "__main__":
    entropy_dir = "F:\\Project_data\\ENTROPY\\Entropy3_fixed_60_wo_plots"
    ss_dir = "F:\\Project_data\\ENTROPY\\SS"  # Update this if the path to SS files is different
    output_dir = "F:\\Project_data\\ENTROPY\\Entropy3_fixed_60_wo_plots"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Debug: List files in entropy directory
    entropy_files = os.listdir(entropy_dir)
    print(f"Files in entropy directory: {entropy_files}")

    for entropy_file_name in entropy_files:
        if entropy_file_name.startswith("entropy_") and entropy_file_name.endswith(".csv"):
            entropy_file = os.path.join(entropy_dir, entropy_file_name)
            protein_name = os.path.basename(entropy_file).replace("entropy_", "").replace(".csv", "")
            ss_file = os.path.join(ss_dir, f"SS_{protein_name}.csv")

            print(f"Processing: {protein_name}")
            print(f"Found entropy file: {entropy_file}")
            print(f"Looking for SS file: {ss_file}")

            if os.path.exists(ss_file):
                plot_entropy_and_ss(entropy_file, ss_file, output_dir, protein_name)
            else:
                print(f"Secondary structure file for {protein_name} not found")
