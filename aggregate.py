import pandas as pd
import os

def aggregate_delta_files(directory, output_directory):
    for sub_dir in os.listdir(directory):
        sub_dir_path = os.path.join(directory, sub_dir)
        
        if os.path.isdir(sub_dir_path):  # Check if it's a directory
            aggregated_data = pd.DataFrame()
            
            for filename in os.listdir(sub_dir_path):
                if filename.endswith('_delta_phi_psi.csv'):
                    file_path = os.path.join(sub_dir_path, filename)
                    
                    # Check if file exists and is not empty
                    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                        try:
                            df = pd.read_csv(file_path, usecols=['RefResidueNumber', 'DeltaPhi', 'DeltaPsi'])
                            df['Filename'] = filename  # Optionally keep track of which file the data came from
                            aggregated_data = pd.concat([aggregated_data, df], ignore_index=True)
                        except pd.errors.EmptyDataError:
                            print(f"Skipping empty or invalid file: {filename}")
                        except ValueError as e:
                            print(f"Error reading {filename}: {e}")
                    else:
                        print(f"File is empty or does not exist: {filename}")
            
            # Generate output filename using the subdirectory name only
            output_file = os.path.join(output_directory, f"{sub_dir}.csv")
            aggregated_data.to_csv(output_file, index=False)
            print("Data aggregation complete and saved to:", output_file)

if __name__ == "__main__":
    # Specify the main directory and the output directory here
    main_dir = "f:\\ISHANI\\Project_data\\delta_filtered\\new"
    output_dir = "f:\\ISHANI\\Project_data\\Entropy\\aggregate\\new"
    aggregate_delta_files(main_dir, output_dir)
