import re
import py3Dmol
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

class NPA:
    def __init__(self, filepath):
        self.filepath = filepath
        self.extract_npa_data()

    def extract_npa_data(self):
        self.npa_data = {}
        with open(self.filepath, 'r') as f:
            lines = f.readlines()

        npa_summary_started = False
        for line in lines:
            if "Summary of Natural Population Analysis:" in line:
                npa_summary_started = True
                continue  

            if npa_summary_started:
                
                if "Atom No    Charge" in line:  
                    continue  

                if "=============" in line or "Total" in line.split() or not line.strip():
                    if self.npa_data:  
                        break
                    else:
                        continue  

                parts = line.split()
                if len(parts) >= 6:
                    try:
                        if re.match(r'^[A-Za-z]+\d+$', parts[0]):
                            atom_label = parts[0]
                            offset = 0
                        elif re.match(r'^[A-Za-z]+$', parts[0]) and re.match(r'^\d+$', parts[1]):
                            atom_label = parts[0] + parts[1]
                            offset = 1
                        else:
                            raise ValueError("Unexpected atom label format")
                        charge = float(parts[offset + 1])   
                        core = float(parts[offset + 2])       
                        valence = float(parts[offset + 3])    
                        rydberg = float(parts[offset + 4])    
                        total = float(parts[offset + 5])
                        spin_density = float(parts[offset + 6]) if len(parts) > offset + 6 else None      

                        self.npa_data[atom_label] = {
                            "Natural Charge": charge,
                            "Core": core,
                            "Valence": valence,
                            "Rydberg": rydberg,
                            "Total": total,
                            "Spin Density": spin_density
                        }
                    except (ValueError, IndexError) as e:
                        print(f"Parse error in line: {line.strip()}")
                continue

        return self.npa_data
    
    # --- print the NPA data ---
    def print_npa_data(self):
        for atom, data in self.npa_data.items():
            print(f"Atom: {atom}")
            for key, value in data.items():
                print(f"  {key}: {value}")

    # --- visualise with py3dmol ---
    def visualise_property(self, xyz_file, property_name = "Natural Charge", gradient = "rainbow", labels=True, stick_size=0.15, sphere_size=0.25):
        g = gradient.lower()
        xyz_data = open(xyz_file).read()
        xyz_lines = xyz_data.splitlines()
        num_atoms_xyz = int(xyz_lines[0])
        print("Number of atoms in XYZ file:", num_atoms_xyz)

        view = py3Dmol.view(width=1000, height=800)
        view.addModel(xyz_data, 'xyz')
        view.setStyle({'stick': {'radius': 0.15}})  # Keep bonds visible

        # --- Color mapping based on NPA value ---
        values = []
        atom_indices_with_value = []
        xyz_elements = [line.split()[0] for line in xyz_data.split('\n')[2:2+num_atoms_xyz]]  # Extract XYZ elements

        for atom_label, data in self.npa_data.items():
            # Parse atom label with updated regex
            match = re.match(r'^([A-Za-z]+)(\d+)$', atom_label, re.IGNORECASE)
            if not match:
                print(f"Warning: Could not parse atom label: {atom_label}")
                continue
            
            element = match.group(1).upper()
            atom_index = int(match.group(2)) - 1  # Convert to 0-based index
            
            # Validate index and element
            if 0 <= atom_index < num_atoms_xyz:
                if xyz_elements[atom_index] != element:
                    print(f"Element mismatch at index {atom_index}: XYZ={xyz_elements[atom_index]}, NPA={element}")
                    continue
                
                try:
                    value = float(data[property_name])
                    values.append(value)
                    atom_indices_with_value.append(atom_index)
                except (ValueError, TypeError):
                    print(f"Invalid value for {atom_label}: {data[property_name]}")
            else:
                print(f"Warning: Atom index {atom_index + 1} from NPA data out of range (XYZ has {num_atoms_xyz} atoms)")

        # if not values:
        #     print("Error: No valid values found. Using default colors.")
        #     view.setStyle({}, {'sphere': {'color': 'grey', 'scale': 0.25}})
        else:
            min_value, max_value = min(values), max(values)
            norm = mcolors.Normalize(vmin=min_value, vmax=max_value)
            cmap = cm.get_cmap(g)

            styled_indices = set()
            for atom_label, data in self.npa_data.items():
                match = re.match(r'^([A-Za-z]+)(\d+)$', atom_label, re.IGNORECASE)
                if not match:
                    continue
                    
                atom_index = int(match.group(2)) - 1
                if 0 <= atom_index < num_atoms_xyz and xyz_elements[atom_index] == match.group(1).upper():
                    try:
                        value = float(data[property_name])
                        hex_color = mcolors.to_hex(cmap(norm(value)))
                        # Maintain both stick and sphere styles
                        view.setStyle({'index': atom_index}, {
                            'stick': {'radius': stick_size},
                            'sphere': {'color': hex_color, 'scale': sphere_size}
                        })
                        styled_indices.add(atom_index)
                    except (ValueError, TypeError):
                        continue

            #  missing atoms (hopefully not needed)
            missing_indices = set(range(num_atoms_xyz)) - styled_indices
            if missing_indices:
                view.setStyle({'index': list(missing_indices)}, {'sphere': {'color': 'grey', 'scale': 0.2}})

            print(f"NPA Value range: {min_value:.3f} to {max_value:.3f}")
            print(f"Found {len(styled_indices)}/{num_atoms_xyz} atoms")
        if labels:
            for atom_label, data in self.npa_data.items():
                match = re.match(r'^([A-Za-z]+)(\d+)$', atom_label, re.IGNORECASE)
                if not match:
                    continue
                element = match.group(1).upper()
                atom_index = int(match.group(2)) - 1  # 0 indexed
                if 0 <= atom_index < num_atoms_xyz and xyz_elements[atom_index] == element:
                    try:
                        value = float(data[property_name])
                    except (ValueError, TypeError):
                        continue
                    parts = xyz_lines[atom_index+2].split()
                    x, y, z = map(float, parts[1:4])

                    hex_color = mcolors.to_hex(cmap(norm(value)))
                    view.addLabel(f"{value:.3f}", {
                        'position': {'x': x, 'y': y, 'z': z},
                        'backgroundColor': hex_color,
                        'backgroundOpacity': 0.5,
                        'fontColor': 'black',
                        'fontSize': 10
                    })

        # --- Color bar for visualization ---
        fig, ax = plt.subplots(figsize=(6, 1))
        fig.subplots_adjust(bottom=0.5)
        cmap = cm.get_cmap(g)
        norm = mcolors.Normalize(vmin=min_value, vmax=max_value)
        cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='horizontal')
        cb.set_label(property_name)
        plt.show()
        view.zoomTo()
        view.show()