import os
import subprocess
import re

# Directory setup
pdb_dir = "/mnt/"  # Your PDB files directory
output_dir = "/mnt/"
reference_pdb_name = ".pdb"  # Name of the reference PDB file

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Path to TM-align executable
tm_align_path = "/usr/local/bin/TMalign"

# Function to parse RMSD from TM-align output
def parse_rmsd(output):
    match = re.search(r"RMSD=\s+(\d+\.\d+)", output)
    if match:
        return float(match.group(1))
    return None

# Run TM-align for each PDB file against the reference and save the output if RMSD <= 3
reference_pdb_path = os.path.join(pdb_dir, reference_pdb_name)

for pdb_file in os.listdir(pdb_dir):
    if pdb_file == reference_pdb_name or not pdb_file.endswith(".pdb"):
        continue
    
    pdb_path = os.path.join(pdb_dir, pdb_file)
    
    cmd = [tm_align_path, pdb_path, reference_pdb_path]
    
    print(f"Running command: {' '.join(cmd)}")  # Debugging statement
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    print(f"TM-align output for {pdb_file}:\n{result.stdout}")  # Debugging statement
    
    rmsd = parse_rmsd(result.stdout)
    
    print(f"RMSD for {pdb_file}: {rmsd}")  # Debugging statement
    
    if rmsd is not None and rmsd <= 3.0:
        output_file_path = os.path.join(output_dir, f"{os.path.splitext(pdb_file)[0]}_vs_{os.path.splitext(reference_pdb_name)[0]}.txt")
        with open(output_file_path, 'w') as output_file:
            output_file.write(result.stdout)
        print(f"Output saved to: {output_file_path}")  # Debugging statement

print("TM-align processing complete.")
