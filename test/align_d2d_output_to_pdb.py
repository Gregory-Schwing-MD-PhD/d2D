import pandas as pd
import MDAnalysis as mda
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from prody import fetchPDB, parsePDB, writePDB

# Your three_to_one dict and function
three_to_one_dict = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}

def three_to_one(resname):
    return three_to_one_dict.get(resname.upper(), 'X')

# Parse SS output file (your code)
def parse_d2d_output(filename):
    records = []
    in_data_section = False
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("#num") and "res" in line:
                in_data_section = True
                continue
            if in_data_section:
                if line.startswith("#") or not line:
                    continue
                parts = line.split()
                if len(parts) < 7:
                    continue
                try:
                    res_num = int(parts[0])
                    residue = parts[1]
                    helix = float(parts[2])
                    beta = float(parts[3])
                    coil = float(parts[4])
                    ppii = float(parts[5])
                    ss = parts[6]
                    records.append({
                        'res_num': res_num,
                        'AA': residue,
                        'Helix': helix,
                        'Beta': beta,
                        'Coil': coil,
                        'PPII': ppii,
                        'SS': ss
                    })
                except ValueError:
                    continue
    return pd.DataFrame.from_records(records)

# Load SS data
df = parse_d2d_output("SS-results.dat")

# Save to CSV without the DataFrame index
df['Helix_normalized'] = (df['Helix'] - df['Helix'].min()) / (df['Helix'].max() - df['Helix'].min())
df.to_csv('d2d.csv', index=False)

seq1 = ''.join(df['AA'].values)

# Fetch PDB with ProDy
pdb_id = "6G6K"
filename = fetchPDB(pdb_id)
structure = parsePDB(filename)
writePDB(filename, structure.protein)  # overwrite with protein only

# Load with MDAnalysis
u = mda.Universe(filename)

# Select chain A
chain_A = u.select_atoms("segid A")

# Function to convert 3-letter to 1-letter safely
def safe_three_to_one(resname):
    return three_to_one(resname.capitalize())

seq2 = ''.join(safe_three_to_one(res.resname) for res in chain_A.residues)

# Align sequences
alignments = pairwise2.align.globalxx(seq1, seq2)
best_alignment = alignments[0]
print(format_alignment(*best_alignment))

aligned_seq1, aligned_seq2, score, begin, end = best_alignment

# Now map helix propensity from df to residues in chain_A using alignment

# Index for df (seq1) and chain_A residues (seq2)
df_index = 0
chain_index = 0

# Create a list to store helix propensities for chain_A residues
helix_propensities = [0.0] * len(chain_A.residues)

for i in range(len(aligned_seq1)):
    aa1 = aligned_seq1[i]
    aa2 = aligned_seq2[i]
    if aa1 != '-' and aa2 != '-':
        # Both aligned residues correspond to actual residues in seq1 and seq2
        # Assign helix propensity from df at df_index to chain_A at chain_index
        #helix_propensities[chain_index] = df.loc[df_index, 'Helix']
        helix_propensities[chain_index] = df.loc[df_index, 'Helix_normalized']
        df_index += 1
        chain_index += 1
    elif aa1 == '-' and aa2 != '-':
        # Gap in seq1, residue in chain_A not aligned -> set helix propensity to 0 or default
        helix_propensities[chain_index] = 0.0
        chain_index += 1
    elif aa1 != '-' and aa2 == '-':
        # Residue in seq1 but gap in chain_A -> skip df_index
        df_index += 1
    else:
        # Both gaps, rare case, skip both
        pass

# Now assign these helix propensities as B-factors to all atoms of each residue in chain_A
for res_idx, res in enumerate(chain_A.residues):
    helix_val = helix_propensities[res_idx]
    # Assign helix_val to B-factor of all atoms in this residue
    res.atoms.tempfactors = helix_val

# Optionally save to new PDB file
with mda.Writer("chainA_with_helix_propensity.pdb", chain_A.n_atoms) as W:
    W.write(chain_A)

print("Helix propensity assigned to B-factors in chain A and saved to 'chainA_with_helix_propensity.pdb'")

import matplotlib.pyplot as plt

# Use beta propensities instead of helix for plotting
helix_protensities = [0.0] * len(chain_A.residues)

df_index = 0
chain_index = 0

for i in range(len(aligned_seq1)):
    aa1 = aligned_seq1[i]
    aa2 = aligned_seq2[i]
    if aa1 != '-' and aa2 != '-':
        helix_protensities[chain_index] = df.loc[df_index, 'Helix_normalized']
        #helix_protensities[chain_index] = df.loc[df_index, 'Helix']
        df_index += 1
        chain_index += 1
    elif aa1 == '-' and aa2 != '-':
        helix_protensities[chain_index] = 0.0
        chain_index += 1
    elif aa1 != '-' and aa2 == '-':
        df_index += 1

import pandas as pd

# Create a DataFrame with residue index and Helix Propensity
df_out = pd.DataFrame({
    'residue_index': range(1, len(helix_protensities) + 1),
    'helix_propensity': helix_protensities
})

# Save to CSV without the DataFrame index
df_out.to_csv('helix_propensity_chainA.csv', index=False)


# Plot
plt.figure(figsize=(12, 4))
plt.plot(range(1, len(helix_protensities) + 1), helix_protensities, label='Helix Propensity', color='red')
plt.xlabel('Residue Index (chain A)')
plt.ylabel('Helix Propensity (B-factor)')
plt.title('Helix Propensity per Residue in Chain A')
plt.legend()
plt.grid(True)
plt.show()

import pymol2
import numpy as np

def align_longest_axis_to_z(cmd, selection="chainA"):
    coords = np.array(cmd.get_coords(selection))
    cov = np.cov(coords.T)
    
    # Eigen decomposition for principal axes
    eigvals, eigvecs = np.linalg.eigh(cov)
    
    # Longest axis = eigenvector with largest eigenvalue
    idx = np.argmax(eigvals)
    longest_axis = eigvecs[:, idx]
    
    # Ensure axis points upward (+z)
    target = np.array([0, 0, 1])
    if np.dot(longest_axis, target) < 0:
        longest_axis = -longest_axis
    
    # Calculate rotation axis and angle
    rot_axis = np.cross(longest_axis, target)
    norm = np.linalg.norm(rot_axis)
    if norm < 1e-6:
        # Already aligned, no rotation needed
        return
    
    rot_axis /= norm
    angle = np.arccos(np.clip(np.dot(longest_axis, target), -1, 1))
    angle_deg = np.degrees(angle)
    
    cmd.rotate(rot_axis.tolist(), angle_deg, selection=selection)

# Run PyMOL session
with pymol2.PyMOL() as pymol:
    cmd = pymol.cmd
    pdb_file = "chainA_with_helix_propensity.pdb"
    
    cmd.load(pdb_file, "chainA")
    cmd.hide("everything", "chainA")
    cmd.show("cartoon", "chainA")
    
    # Center and orient first
    cmd.orient("chainA")
    
    # Align longest principal axis with +Z
    align_longest_axis_to_z(cmd, "chainA")
    
    # Color by B-factor (Helix Propensity)
    cmd.spectrum("b", "blue_white_red", "chainA", minimum=0, maximum=1)
    
    cmd.zoom("chainA")
    cmd.rotate("z", -90, "chainA")  # rotate 90 degrees about x axis
    # Assuming Helix Propensity values range from 0 (min) to 1 (max)
    cmd.spectrum("b", "blue_white_red", "chainA", minimum=0, maximum=1)
    # Save image
    cmd.png("chainA_with_helix_propensity.png", width=800, height=600, dpi=300, ray=1)
