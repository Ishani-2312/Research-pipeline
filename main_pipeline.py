import os
import extract_chainA  # Code 1
import calculate_phi_psi  # Code 2
import run_tm_align  # Code 3
import residue_mapping  # Code 4
import delta_phi_psi  # Code 5
import entropy_calculation  # Code 6
import run_dssp  # Code 7
import b_factor_extraction  # Code 8
import rmsd_calculation  # Code 9
import asa_extraction  # Code 10

def main_pipeline(input_dir, output_dir):
    print("Starting the pipeline...")

    # Step 1: Extract Chain A from PDB files
    chainA_output_dir = os.path.join(output_dir, "pdb_chainA")
    extract_chainA.process_structures_and_save_angles(input_dir, chainA_output_dir)
    
    # Step 2: Calculate phi/psi angles
    dihedral_output_dir = os.path.join(output_dir, "dihedral")
    calculate_phi_psi.process_structures_and_save_angles(chainA_output_dir, dihedral_output_dir)
    
    # Step 3: Run TM-align
    tm_output_dir = os.path.join(output_dir, "TM_output")
    run_tm_align.run_tm_align(chainA_output_dir, tm_output_dir)
    
    # Step 4: Compare residue mappings
    residue_mapping_dir = os.path.join(output_dir, "residue_mapping")
    residue_mapping.process_residue_mappings(tm_output_dir, chainA_output_dir, residue_mapping_dir)
    
    # Step 5: Calculate Delta Phi and Psi
    delta_output_dir = os.path.join(output_dir, "delta")
    delta_phi_psi.process_alignment_pairs(residue_mapping_dir, dihedral_output_dir, delta_output_dir)
    
    # Step 6: Calculate entropy
    entropy_output_dir = os.path.join(output_dir, "entropy")
    entropy_calculation.calculate_entropy(delta_output_dir, entropy_output_dir)
    
    # Step 7: Run DSSP for secondary structure analysis
    dssp_output_dir = os.path.join(output_dir, "dssp_analysis")
    run_dssp.process_pdb_files(chainA_output_dir, dssp_output_dir)
    
    # Step 8: Extract B-factors
    b_factor_output_dir = os.path.join(output_dir, "b_factors")
    b_factor_extraction.process_pdb_files(chainA_output_dir, b_factor_output_dir)
    
    # Step 9: Calculate RMSD
    rmsd_output_dir = os.path.join(output_dir, "rmsd")
    rmsd_calculation.process_rmsd_files(chainA_output_dir, rmsd_output_dir)
    
    # Step 10: Extract ASA
    asa_output_dir = os.path.join(output_dir, "asa")
    asa_extraction.process_dssp_files(dssp_output_dir, asa_output_dir)

    print("Pipeline completed!")

if __name__ == "__main__":
    input_directory = "your/input/directory/path"  # Replace with your input directory
    output_directory = "your/output/directory/path"  # Replace with your output directory
    main_pipeline(input_directory, output_directory)
