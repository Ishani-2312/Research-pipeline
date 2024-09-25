import os
import pandas as pd
import matplotlib.pyplot as plt

def extract_percentages(csv_file):
    """
    Extracts secondary structure percentages from a given CSV file.
    Returns a DataFrame with SecondaryStructure and Percentage columns.
    """
    try:
        with open(csv_file, 'r') as f:
            lines = f.readlines()
            # Find the start of the "Secondary Structure Percentages" section
            start_line = -1
            for i, line in enumerate(lines):
                if "Secondary Structure Percentages" in line:
                    start_line = i + 1
                    break
            if start_line == -1:
                print(f"No percentages section found in {csv_file}")
                return pd.DataFrame()

        # Extract the relevant section directly using pandas
        ss_percentages = pd.read_csv(csv_file, skiprows=start_line, header=None)
        ss_percentages.columns = ['SecondaryStructure', 'Percentage']

        # Convert 'Percentage' to numeric, coercing errors to NaN and then filling NaNs with zero
        ss_percentages['Percentage'] = pd.to_numeric(ss_percentages['Percentage'], errors='coerce').fillna(0)

        # Filter out invalid rows
        ss_percentages = ss_percentages[ss_percentages['SecondaryStructure'].isin(['Helix', 'Beta Strand', 'Turn', 'Coil'])]

        # Add missing secondary structure types with zero percentage if not present
        ss_types = ['Helix', 'Beta Strand', 'Turn', 'Coil']
        existing_ss = ss_percentages['SecondaryStructure'].values
        missing_ss = [{'SecondaryStructure': ss_type, 'Percentage': 0} for ss_type in ss_types if ss_type not in existing_ss]
        if missing_ss:
            ss_percentages = pd.concat([ss_percentages, pd.DataFrame(missing_ss)], ignore_index=True)

        # Debug: Print the extracted DataFrame
        print(f"Extracted from {csv_file}:")
        print(ss_percentages)
        return ss_percentages
    except pd.errors.EmptyDataError:
        print(f"Secondary structure percentages not found in {csv_file}")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error processing {csv_file}: {e}")
        return pd.DataFrame()

def aggregate_percentages(input_directory):
    """
    Aggregates secondary structure percentages from all CSV files in the input directory.
    Ensures that all secondary structure types are represented for each protein family.
    """
    aggregated_data = []
    ss_types = ['Helix', 'Beta Strand', 'Turn', 'Coil']

    for filename in os.listdir(input_directory):
        if filename.startswith("SS_") and filename.endswith(".csv"):
            csv_file_path = os.path.join(input_directory, filename)
            protein_name = os.path.splitext(filename)[0][3:]  # Strip 'SS_' prefix
            
            ss_percentages = extract_percentages(csv_file_path)
            if not ss_percentages.empty:
                ss_percentages['ProteinFamily'] = protein_name
                aggregated_data.append(ss_percentages)
    
    if aggregated_data:
        # Concatenate all the data
        aggregated_df = pd.concat(aggregated_data, ignore_index=True)
        
        # Ensure each protein family has all SS types
        complete_data = []
        for protein_family in aggregated_df['ProteinFamily'].unique():
            family_data = aggregated_df[aggregated_df['ProteinFamily'] == protein_family]
            for ss_type in ss_types:
                row = family_data[family_data['SecondaryStructure'] == ss_type]
                if row.empty:
                    complete_data.append({'ProteinFamily': protein_family, 'SecondaryStructure': ss_type, 'Percentage': 0})
                else:
                    complete_data.append(row.iloc[0].to_dict())
        
        complete_df = pd.DataFrame(complete_data)
        
        # Debug: Print the aggregated DataFrame
        print("Aggregated DataFrame:")
        print(complete_df)
        return complete_df
    else:
        return pd.DataFrame()

def plot_percentages(data, output_file):
    """
    Plots the secondary structure percentages for each protein family.
    """
    if data.empty:
        print("No data available to plot.")
        return

    # Pivot the data for plotting
    pivot_df = data.pivot(index='ProteinFamily', columns='SecondaryStructure', values='Percentage').fillna(0)

    # Define colors for each secondary structure type
    colors = {
        'Helix': 'red',
        'Beta Strand': 'darkblue',  # Darker blue
        'Turn': 'yellow',
        'Coil': 'green'
    }

    # Plot with specified colors and increased bar width
    ax = pivot_df.plot(kind='bar', figsize=(15, 10), color=[colors[col] for col in pivot_df.columns], width=0.85)

    plt.xlabel('Protein Family', fontsize=18)  # Increased x-axis label size
    plt.ylabel('Percentage', fontsize=18)      # Increased y-axis label size
    plt.title('Secondary Structure Percentages by Protein Family', fontsize=20)
    plt.legend(title='Secondary Structure', fontsize=14)  # Increase legend font size
    plt.xticks(fontsize=16)  # Increased x-axis ticks size, removed rotation
    plt.yticks(fontsize=16)  # Increased y-axis ticks size
    plt.tight_layout()
    
    plt.savefig(output_file)
    plt.show()

if __name__ == "__main__":
    input_directory = "C:\\Users\\msc2\\Desktop\\Project_data\\ENTROPY\\SS"  # Update with your CSV files directory
    output_image_file = "C:\\Users\\msc2\\Desktop\\Project_data\\ENTROPY\\SS_percentages_bar_chart.png"  # Update with desired output image file path
    
    data = aggregate_percentages(input_directory)
    plot_percentages(data, output_image_file)
