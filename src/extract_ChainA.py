import os
import gzip
import shutil
import pymol
from pymol import cmd

# Specify the input and output directories
input_directory = ""
output_directory = ""
os.makedirs(output_directory, exist_ok=True)

# Decompress .gz files
for filename in os.listdir(input_directory):
    if filename.endswith(".gz"):
        filepath = os.path.join(input_directory, filename)
        decompressed_filepath = os.path.join(input_directory, filename[:-3])  # Remove the .gz extension
        with gzip.open(filepath, 'rb') as f_in:
            with open(decompressed_filepath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Decompressed: {filename} -> {decompressed_filepath}")

# Initialize PyMOL
pymol.finish_launching()

# Load decompressed PDB files into PyMOL
for filename in os.listdir(input_directory):
    if filename.endswith(".pdb"):
        filepath = os.path.join(input_directory, filename)
        cmd.load(filepath)

# List all loaded objects in PyMOL
loaded_objects = cmd.get_names("objects")

for object_name in loaded_objects:
    try:
        # Get the list of chains in the current object
        chains = cmd.get_chains(object_name)
        if chains:
            # Select the first chain
            first_chain = chains[0]
            cmd.select('first_chain', f'{object_name} and chain {first_chain}')
            new_filename = f'{object_name}_Chain{first_chain}.pdb'
            new_filepath = os.path.join(output_directory, new_filename)
            cmd.save(new_filepath, 'first_chain')
            print(f"Processed: {object_name} -> {new_filename}")
        else:
            print(f"No chains found in {object_name}")
        
        # Deselect the selection for the next iteration
        cmd.deselect()
    except pymol.CmdException as e:
        print(f"Error processing {object_name}: {e}")

# Display "done" message
print("Done")
