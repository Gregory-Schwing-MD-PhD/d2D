import matplotlib.pyplot as plt
import sys
import os

def parse_d2d_output(filename):
    data = {
        'res_num': [],
        'residue': [],
        'Helix': [],
        'Beta': [],
        'Coil': [],
        'PPII': [],
        'SS': []
    }

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

                    data['res_num'].append(res_num)
                    data['residue'].append(residue)
                    data['Helix'].append(helix)
                    data['Beta'].append(beta)
                    data['Coil'].append(coil)
                    data['PPII'].append(ppii)
                    data['SS'].append(ss)
                except ValueError:
                    continue  # Skip lines with parsing issues

    return data

def plot_secondary_structure(data, filename):
    plt.figure(figsize=(12, 6))
    plt.plot(data['res_num'], data['Helix'], label='Helix', color='red')
    plt.plot(data['res_num'], data['Beta'], label='Beta', color='blue')
    plt.plot(data['res_num'], data['Coil'], label='Coil', color='green')
    plt.plot(data['res_num'], data['PPII'], label='PPII', color='purple')

    plt.xlabel('Residue Number')
    plt.ylabel('Population')
    plt.title(f'Secondary Structure Populations\n{os.path.basename(filename)}')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_d2d.py <d2d_output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    data = parse_d2d_output(input_file)
    plot_secondary_structure(data, input_file)

