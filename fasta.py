from Bio import PDB
import os

def extract_sequence(pdb_file, chain_id='A'):
    """Extracts the sequence for a specified chain from a PDB file."""
    parser = PDB.PDBParser(QUIET=True)  # QUIET=True to avoid unnecessary warnings
    structure = parser.get_structure('PDB_structure', pdb_file)
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                seq = ''
                for residue in chain:
                    if PDB.is_aa(residue, standard=True):
                        seq += PDB.Polypeptide.three_to_one(residue.resname).upper()
                return seq
    return ''

def append_to_fasta(output_file, header, sequence):
    """Appends a sequence to a FASTA file."""
    with open(output_file, 'a') as f:  # Note the 'a' mode for appending
        f.write(f">{header}\n{sequence}\n")

def process_pdb_directory_to_single_fasta(pdb_dir, output_file, chain_id='A'):
    """Processes all PDB files in a directory and appends sequences to a single FASTA file."""
    with open(output_file, 'w'):  # Create or clear the file
        pass
    for filename in os.listdir(pdb_dir):
        if filename.endswith(".pdb"):
            pdb_file = os.path.join(pdb_dir, filename)
            sequence = extract_sequence(pdb_file, chain_id)
            if sequence:
                header = filename[:-4]  # Remove the .pdb extension for the header
                append_to_fasta(output_file, header, sequence)

# Example usage
pdb_directory = "C:\\Users\\msc2\\Desktop\\NEWPIPE\\carbonic_anhydrase\\pdb_chainA"
output_fasta = "C:\\Users\\msc2\\Desktop\\NEWPIPE\\carbonic_anhydrase\\CA.fasta"
process_pdb_directory_to_single_fasta(pdb_directory, output_fasta, 'A')


