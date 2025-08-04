import pandas as pd
import sys

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
                    continue  # Skip incomplete lines

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
                        'residue': residue,
                        'Helix': helix,
                        'Beta': beta,
                        'Coil': coil,
                        'PPII': ppii,
                        'SS': ss
                    })
                except ValueError:
                    continue  # Skip lines with parsing issues

    return pd.DataFrame.from_records(records)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python join_d2d_with_shifty.py <d2d_output_file> <shifty_pickle_file>")
        sys.exit(1)

    d2d_file = sys.argv[1]
    shifty_pickle_file = sys.argv[2]

    d2d_df = parse_d2d_output(d2d_file)
    shifty_df = pd.read_pickle(shifty_pickle_file)

    # Rename '#NUM' to 'res_num' to prepare for merge
    shifty_df = shifty_df.rename(columns={'#NUM': 'res_num'})

    # Merge on 'res_num'
    merged_df = pd.merge(shifty_df, d2d_df, on='res_num', how='inner')

    print(merged_df)
