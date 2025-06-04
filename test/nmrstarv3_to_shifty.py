import sys
import nmrstarlib
import pandas as pd

# Mapping from 3-letter to 1-letter amino acid codes
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XLE': 'J',
    'UNK': 'X'
}

def load_nmrstar_to_df(nmrstar_path):
    starfile = next(nmrstarlib.read_files(nmrstar_path))
    chem_shift_list = starfile['save_assigned_chem_shift_list_1']
    loop_rows = chem_shift_list['loop_1'][1]  # List of OrderedDicts
    df = pd.DataFrame(loop_rows)
    return df

def convert_to_shifty_format(df):
    # Ensure chemical shift values are numeric
    df['Atom_chem_shift.Val'] = pd.to_numeric(df['Atom_chem_shift.Val'], errors='coerce')
    
    # Convert Comp_index_ID to integers for proper grouping
    df['Atom_chem_shift.Comp_index_ID'] = df['Atom_chem_shift.Comp_index_ID'].astype(int)

    # Group without sorting to preserve input order
    grouped = df.groupby('Atom_chem_shift.Comp_index_ID', sort=False)
    
    shifty_rows = []

    for i, (res_id, group) in enumerate(grouped, start=1):
        three_letter = group["Atom_chem_shift.Comp_ID"].iloc[0].upper()
        one_letter = three_to_one.get(three_letter, 'X')  # Default to 'X' if unknown

        row = {
            "#NUM": i,
            "AA": one_letter,
            "HA": 0.0,
            "CA": 0.0,
            "CB": 0.0,
            "CO": 0.0,
            "N": 0.0,
            "HN": 0.0
        }

        for _, atom in group.iterrows():
            atom_name = atom["Atom_chem_shift.Atom_ID"].strip().upper()
            value = atom["Atom_chem_shift.Val"]

            if atom_name == "CA":
                row["CA"] = value
            elif atom_name == "CB":
                row["CB"] = value
            elif atom_name == "C":
                row["CO"] = value
            elif atom_name == "N":
                row["N"] = value
            elif atom_name in {"H", "HN"}:
                row["HN"] = value

        shifty_rows.append(row)

    shifty_df = pd.DataFrame(shifty_rows)
    return shifty_df

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_filename.str> <output_filename>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    df = load_nmrstar_to_df(input_file)
    shifty_df = convert_to_shifty_format(df)

    shifty_df.to_csv(output_file, index=False, sep="\t")
    print(f"Shifty format data exported to: {output_file} (tab-delimited)")

if __name__ == "__main__":
    main()

