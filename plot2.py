import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sklearn.preprocessing import MinMaxScaler

def scale_data(df, columns):
    """Scale the specified columns of a DataFrame using Min-Max scaling."""
    scaler = MinMaxScaler()
    df[columns] = scaler.fit_transform(df[columns])
    return df

def plot_combined_metrics(entropy_file, ss_file, bfactor_file, rmsd_file, sasa_file, output_dir, protein_name):
    try:
        # Load the entropy data
        entropy_data = pd.read_csv(entropy_file)
        print(f"Entropy Data for {protein_name}:")
        print(entropy_data.head())

        # Load the secondary structure data, ignoring the last row if it contains summary information
        ss_data = pd.read_csv(ss_file)
        if isinstance(ss_data.iloc[-1, 0], str) and ss_data.iloc[-1, 0].startswith('Total'):
            ss_data = ss_data.iloc[:-1]
        print(f"SS Data for {protein_name}:")
        print(ss_data.head())

        # Load B-factor, RMSD, and SASA data
        bfactor_data = pd.read_csv(bfactor_file)
        print(f"B-factor Data for {protein_name}:")
        print(bfactor_data.head())
        
        rmsd_data = pd.read_csv(rmsd_file)
        print(f"RMSD Data for {protein_name}:")
        print(rmsd_data.head())
        
        sasa_data = pd.read_csv(sasa_file)
        print(f"SASA Data for {protein_name}:")
        print(sasa_data.head())

        # Rename columns for consistency
        entropy_data.rename(columns={'RefResidueNumber': 'ResidueNumber'}, inplace=True)
        rmsd_data.rename(columns={'Residue': 'ResidueNumber'}, inplace=True)
        sasa_data.rename(columns={'Residue_Number': 'ResidueNumber', 'ASA': 'SASA'}, inplace=True)

        # Ensure the columns to be merged have the same data type
        entropy_data['ResidueNumber'] = entropy_data['ResidueNumber'].astype(str)
        ss_data['ResidueNumber'] = ss_data['ResidueNumber'].astype(str)
        bfactor_data['ResidueNumber'] = bfactor_data['ResidueNumber'].astype(str)
        rmsd_data['ResidueNumber'] = rmsd_data['ResidueNumber'].astype(str)
        sasa_data['ResidueNumber'] = sasa_data['ResidueNumber'].astype(str)

        # Merge all datasets on the residue number
        merged_data = pd.merge(entropy_data, ss_data, on='ResidueNumber')
        merged_data = pd.merge(merged_data, bfactor_data, on='ResidueNumber')
        merged_data = pd.merge(merged_data, rmsd_data, on='ResidueNumber')
        merged_data = pd.merge(merged_data, sasa_data, on='ResidueNumber')
        
        print(f"Merged Data for {protein_name}:")
        print(merged_data.head())

        # Convert merged residue numbers to numeric for plotting
        merged_data['ResidueNumber'] = pd.to_numeric(merged_data['ResidueNumber'], errors='coerce')

        # Scale the B-factor, RMSD, and SASA columns
        merged_data = scale_data(merged_data, ['AverageBFactor', 'Mean_RMSD', 'SASA'])

        # Plotting
        fig, (ax1, ax_ss, ax2) = plt.subplots(3, 1, figsize=(28, 12), sharex=True, gridspec_kw={'height_ratios': [5, 0.6, 5]})

        # Plot bar plots for entropy
        bar_width = 0.6  # Increased bar width
        ax1.bar(merged_data['ResidueNumber'] - bar_width/2, merged_data['PhiEntropy'], bar_width, label='Phi Entropy', color='#1F77B4')
        ax1.bar(merged_data['ResidueNumber'] + bar_width/2, merged_data['PsiEntropy'], bar_width, label='Psi Entropy', color='#D62728')

        # Add secondary structure rectangles for helices with reduced height and width
        for i, row in merged_data.iterrows():
            if row['SecondaryStructure'] == 'Helix':
                ax_ss.add_patch(Rectangle((row['ResidueNumber'] - 0.25, 0.2), 0.5, 0.6, color='#D62728'))  # Adjusted width and height
            elif row['SecondaryStructure'] == 'Beta Strand':
                ax_ss.plot(row['ResidueNumber'], 0.5, marker='>', color='#1F77B4', markersize=10, linestyle='None')
            elif row['SecondaryStructure'] in ['Coil', 'Turn']:
                ax_ss.plot(row['ResidueNumber'], 0.5, marker='_', color='#2CA02C', markersize=10, linestyle='None')

        # Remove y-axis and ticks for secondary structure plot
        ax_ss.set_yticks([])
        ax_ss.set_yticklabels([])
        ax_ss.set_ylim(0, 1)
        ax_ss.set_ylabel('SS', fontsize=18)

        # Second subplot for additional metrics
        ax2.plot(merged_data['ResidueNumber'], merged_data['AverageBFactor'], label='B-Factor', color='#4B0082', linewidth=2, linestyle='--')
        ax2.plot(merged_data['ResidueNumber'], merged_data['Mean_RMSD'], label='Mean RMSD', color='#556B2F', linewidth=2, linestyle='-.')
        ax2.plot(merged_data['ResidueNumber'], merged_data['SASA'], label='SASA', color='#FF1493', linewidth=2, linestyle=':')

        # Set y-axis limits for entropy
        ax1.set_ylim(0, 5)

        # Set labels
        ax2.set_xlabel('Residue Number', fontsize=20)
        ax1.set_ylabel('Entropy', fontsize=20)
        ax2.set_ylabel('B-Factor, RMSD, SASA', fontsize=20)

        # Increasing font sizes for ticks
        ax1.tick_params(axis='both', which='major', labelsize=20)
        ax2.tick_params(axis='both', which='major', labelsize=20)

        # Set x-ticks at regular intervals to avoid clutter
        interval = max(1, len(merged_data) // 20)
        xticks = merged_data['ResidueNumber'][::interval]
        ax2.set_xticks(xticks)
        ax2.set_xticklabels(xticks, rotation=45, ha='right')

        # Create a secondary x-axis below the main plot to display residue numbers
        ax_resnum = ax1.secondary_xaxis('bottom')  # Position closer to the main x-axis
        ax_resnum.set_xlim(ax1.get_xlim())
        ax_resnum.set_xticks(xticks)
        ax_resnum.set_xticklabels(xticks)
        ax_resnum.tick_params(axis='x', labelsize=18, rotation=45)  # Rotate by 45 degrees and set font size
 
        # Add grid lines
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)

        plt.subplots_adjust(top=0.95, bottom=0.1, hspace=0.3)  # Adjust the top and bottom margins and the height spacing

        # Save the plot to a file with high quality
        output_file = os.path.join(output_dir, f"combined_plot_{protein_name}.png")
        plt.savefig(output_file, dpi=600)

        # Show the plot (optional)
        plt.show()
        print(f"Combined plot generated and saved for {protein_name}")

    except Exception as e:
        print(f"An error occurred while processing {protein_name}: {e}")

if __name__ == "__main__":
    entropy_dir = "F:\\Project_data\\ENTROPY\\Entropy3_fixed_60_wo_plots"
    ss_dir = "F:\\Project_data\\ENTROPY\\SS"
    bfactor_dir = "F:\\Project_data\\ENTROPY\\b-factor"
    rmsd_dir = "F:\\Project_data\\ENTROPY\\rmsd"
    sasa_dir = "F:\\Project_data\\ASA_wo_norm"
    output_dir = "F:\\Project_data\\ENTROPY\\combine_plots"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Debug: List files in entropy directory
    entropy_files = os.listdir(entropy_dir)
    ss_files = os.listdir(ss_dir)
    bfactor_files = os.listdir(bfactor_dir)
    rmsd_files = os.listdir(rmsd_dir)
    sasa_files = os.listdir(sasa_dir)

    print(f"Files in entropy directory: {entropy_files}")
    print(f"Files in SS directory: {ss_files}")
    print(f"Files in B-factor directory: {bfactor_files}")
    print(f"Files in RMSD directory: {rmsd_files}")
    print(f"Files in SASA directory: {sasa_files}")

    for entropy_file_name in entropy_files:
        if entropy_file_name.startswith("entropy_") and entropy_file_name.endswith(".csv"):
            entropy_file = os.path.join(entropy_dir, entropy_file_name)
            protein_name = os.path.basename(entropy_file).replace("entropy_", "").replace(".csv", "")
            ss_file = os.path.join(ss_dir, f"SS_{protein_name}.csv")
            bfactor_file = os.path.join(bfactor_dir, f"{protein_name}.csv")
            rmsd_file = os.path.join(rmsd_dir, f"results_mean_RMSD_{protein_name}.csv")
            sasa_file = os.path.join(sasa_dir, f"{protein_name}_sasa.csv")

            print(f"Processing: {protein_name}")
            print(f"Found entropy file: {entropy_file}")
            print(f"Looking for SS file: {ss_file}")
            print(f"Looking for B-factor file: {bfactor_file}")
            print(f"Looking for RMSD file: {rmsd_file}")
            print(f"Looking for SASA file: {sasa_file}")

            if os.path.exists(ss_file) and os.path.exists(bfactor_file) and os.path.exists(rmsd_file) and os.path.exists(sasa_file):
                plot_combined_metrics(entropy_file, ss_file, bfactor_file, rmsd_file, sasa_file, output_dir, protein_name)
            else:
                print(f"Required files for {protein_name} not found")
