import sys
import pandas as pd
import pynmrstar
from io import StringIO

# Mapping from 3-letter to 1-letter amino acid codes
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XLE': 'J',
    'UNK': 'X'
}

def load_nmrstar_to_df(entry_id):
    print(f"📥 Downloading BMRB entry {entry_id}...")
    entry = pynmrstar.Entry.from_database(entry_id)
    csv = entry['assigned_chem_shift_list_1']['_Atom_chem_shift'].get_data_as_csv()
    df = pd.read_csv(StringIO(csv))
    return df

def convert_to_shifty_format(df):
    # Clean column names
    df.columns = [col.strip() for col in df.columns]

    # Rename for convenience
    df = df.rename(columns={
        '_Atom_chem_shift.Comp_index_ID': 'ResID',
        '_Atom_chem_shift.Comp_ID': 'ResName',
        '_Atom_chem_shift.Atom_ID': 'Atom',
        '_Atom_chem_shift.Val': 'Value'
    })

    df = df[['ResID', 'ResName', 'Atom', 'Value']].copy()

    # Normalize atom names
    df['Atom'] = df['Atom'].str.strip().str.upper()
    df['Atom'] = df['Atom'].replace({'H': 'HN'})

    # Pivot into wide format
    pivot = df.pivot_table(
        index=['ResID', 'ResName'],
        columns='Atom',
        values='Value',
        aggfunc='first'
    ).reset_index()

    # Rename atoms to Shifty column names
    pivot = pivot.rename(columns={'C': 'CO', 'HN': 'HN'})

    # Ensure required atoms are present
    for col in ['CA', 'CB', 'CO', 'N', 'HN']:
        if col not in pivot.columns:
            pivot[col] = 0.0

    # Map 3-letter to 1-letter amino acids
    pivot['AA'] = pivot['ResName'].str.upper().map(three_to_one).fillna('X')

    # Rename and sort
    pivot = pivot.rename(columns={'ResID': '#NUM'})
    pivot['#NUM'] = pivot['#NUM'].astype(int)
    pivot = pivot.sort_values('#NUM')

    # Add HA if missing
    if 'HA' not in pivot.columns:
        pivot['HA'] = 0.0

    # Final order
    pivot = pivot[['#NUM', 'AA', 'HA', 'CA', 'CB', 'CO', 'N', 'HN']]

    # Fill any remaining missing values with 0.0
    pivot = pivot.fillna(0.0)

    return pivot

def main():
    if len(sys.argv) != 3:
        print("Usage: python shifty_exporter.py <BMRB_entry_ID> <output_file.tsv>")
        sys.exit(1)

    entry_id = sys.argv[1]
    output_file = sys.argv[2]

    try:
        df = load_nmrstar_to_df(entry_id)
        shifty_df = convert_to_shifty_format(df)
        shifty_df.to_csv(output_file, index=False, sep="\t")
        print(f"✅ Shifty format data exported to: {output_file}")
        shifty_df.to_pickle(output_file+".pkl")
        print(f"✅ Data frame exported to: {output_file}.pkl")
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

