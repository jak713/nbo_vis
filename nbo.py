from matplotlib import cm, colors
import re
import py3Dmol
import matplotlib.pyplot as plt

# SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS

#      Threshold for printing:   0.50 kcal/mol
#     (Intermolecular threshold: 0.05 kcal/mol)
#                                                           E(2) E(NL)-E(L) F(L,NL)
#       Donor (L) NBO              Acceptor (NL) NBO      kcal/mol   a.u.    a.u.
#  ===============================================================================
                
#  within unit  1
#   153. LP ( 1) N  1           555. BD*( 1) C 30- C 52     11.15    1.09   0.099
#   153. LP ( 1) N  1           841. RY ( 9) N  1            0.70    4.88   0.052
#   153. LP ( 1) N  1          1398. RY ( 1) C 30           13.80    1.57   0.132
#   160. LP ( 1) C116           649. BD*( 2) C 96- C103    146.42    0.20   0.154
#   160. LP ( 1) C116           651. BD*( 2) C 97- C104    139.08    0.21   0.153
#   160. LP ( 1) C116           718. BD*( 1) C133- C153      6.64    0.76   0.063
#   160. LP ( 1) C116           719. BD*( 1) C133- C154      6.02    0.76   0.060
#   160. LP ( 1) C116           720. BD*( 2) C133- C154     35.66    0.20   0.076


### LV = lone valence orbital, very low occupancy 
# Donor (donates electron density) = LP or BD (L=Lewis) if lone pair then the donor is one atom if bonding orbital then the donor is two atoms connected by a bond
#Acceptor (receives electron density hence being stabilised) = BD* or RY (NL=Non-Lewis) if rydberg orbital then the acceptor is one atom if antibodning orbital then the acceptor is two atoms connected by a bond
# the example data to be extracted is as follows:
#  153. LP ( 1) N  1           555. BD*( 1) C 30- C 52     11.15    1.09   0.099
# 153 = specific NBO index for the specific orbital (important for the NBO analysis)
# LP = lone pair donor
# ( 1) = the first lone pair on the atom
# N = nitrogen atom
# 1 = atom number
# 555 = specific NBO index for the specific orbital (important for the NBO analysis)
# BD* = anti-bonding orbital acceptor
# ( 1) = the first anti-bonding orbital between the two atoms
# C 30- C 52 = the two atoms connected by the bond
# 11.15 = E(2) which is the stabilisation energy in kcal/mol
# 1.09 = the energy difference between the donor and the acceptor orbitals in a.u.
# 0.099 = the fock matrix element between the donor and the acceptor orbitals in a.u. (quantifies interaction strength)

# need to design the right kind of regex to extract the data from the lines and the right kind of dictionary to store the data
class NBO_SOP:
    def __init__(self, filepath):
        self.filepath = filepath
        self.extract_nbo_data()

    def extract_nbo_data(self):
        self.nbo_data = []
        with open(self.filepath, 'r') as f:
            lines = f.readlines()

        nbo_analysis_started = False
        for line in lines:
            if "SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS" in line:
                nbo_analysis_started = True
                continue  
    
            if nbo_analysis_started:
                if "Donor (L) NBO" in line:  
                    continue  
                if line.strip() == "":
                    continue
                if "from" in line:
                    continue
                if "within unit" in line:
                    continue
                if "threshold" in line or "Threshold" in line:
                    continue
                if "E(2) E(NL)-E(L) F(L,NL)" in line:
                    continue
                if "=" in line:
                    continue
                if "NATURAL BOND ORBITALS (Summary):" in line:
                    break


                pattern = re.compile(r"""
                    ^\s*
                    (\d+)\.\s+                          # Donor Index (1)
                    (LP|BD)\s*\(\s*(\d+)\s*\)\s+       # Donor Type (2), Donor Orb No (3)
                    (.+?)\s{3,}                         # Donor Atom String (4) - Capture until 2+ spaces (separator)
                    (\d+)\.\s+                          # Acceptor Index (5)
                    (BD\*|RY|LV)\s*\(\s*(\d+)\s*\)\s+     # Acceptor Type (6), Acceptor Orb No (7)
                    (.+?)\s+                            # Acceptor Atom String (8) - Capture until space before numbers
                    ([\d.-]+)\s+                        # E(2) (9)
                    ([\d.-]+)\s+                        # E Diff (10)
                    ([\d.-]+)                           # Fock Elem (11)
                    \s*$
                """, re.VERBOSE)

                match = pattern.match(line)

                if match:
                    donor_index = match.group(1)
                    donor_type = match.group(2)
                    donor_orb_no = match.group(3)
                    donor_atom_string = match.group(4).strip()
                    acceptor_index = match.group(5)
                    acceptor_type = match.group(6)
                    acceptor_orb_no = match.group(7)
                    acceptor_atom_string = match.group(8).strip()
                    e2 = float(match.group(9))
                    e_diff = float(match.group(10))
                    fock_elem = float(match.group(11))
                    # Split the atom strings into individual atoms
                    donor_atoms = re.split(r'\s*-\s*', donor_atom_string)
                    acceptor_atoms = re.split(r'\s*-\s*', acceptor_atom_string)
                    # Create a dictionary for the NBO data
                else:
                    print(f"Parse error in line: {line.strip()}")
                    continue

                nbo_entry = {
                    "Donor Index": donor_index,
                    "Donor Type": donor_type,
                    "Donor Orb No": donor_orb_no,
                    "Donor Atoms": donor_atoms,
                    "Acceptor Index": acceptor_index,
                    "Acceptor Type": acceptor_type,
                    "Acceptor Orb No": acceptor_orb_no,
                    "Acceptor Atoms": acceptor_atoms,
                    "E(2)": e2,
                    "E Diff": e_diff,
                    "Fock Elem": fock_elem
                }
                self.nbo_data.append(nbo_entry)
        return self.nbo_data
    
    def print_nbo_data(self):
        """Print the NBO data in a formatted table with better alignment"""
        header = (
            f"{'Donor Index':<12} {'Donor Type':<10} {'Donor Orb No':<12} "
            f"{'Donor Atoms':<25} {'Acceptor Index':<15} {'Acceptor Type':<12} "
            f"{'Acceptor Orb No':<15} {'Acceptor Atoms':<25} {'E(2)':>8} {'E Diff':>8} {'Fock Elem':>10}"
        )
        print(header)
        print("=" * len(header))
        for entry in self.nbo_data:
            donor_atoms = ", ".join(entry["Donor Atoms"])
            acceptor_atoms = ", ".join(entry["Acceptor Atoms"])
            row = (
                f"{entry['Donor Index']:<12} {entry['Donor Type']:<10} {entry['Donor Orb No']:<12} "
                f"{donor_atoms:<25} {entry['Acceptor Index']:<15} {entry['Acceptor Type']:<12} "
                f"{entry['Acceptor Orb No']:<15} {acceptor_atoms:<25} {entry['E(2)']:>8.2f} "
                f"{entry['E Diff']:>8.2f} {entry['Fock Elem']:>10.2f}"
            )
            print(row)
        print("=" * len(header))
        print(f"Total number of NBO interactions: {len(self.nbo_data)}")

    def print_loneToAnti(self):
        """Print interactions of donor LP and acceptor BD*"""
        header = (
            f"{'Donor Index':<12} {'Donor Type':<10} {'Donor Orb No':<12} "
            f"{'Donor Atoms':<25} {'Acceptor Index':<15} {'Acceptor Type':<12} "
            f"{'Acceptor Orb No':<15} {'Acceptor Atoms':<25} {'E(2)':>8} {'E Diff':>8} {'Fock Elem':>10}"
        )
        print(header)
        print("=" * len(header))
        counter = 0
        for entry in self.nbo_data:
            if entry["Donor Type"] == "LP" and entry["Acceptor Type"] == "BD*":
                donor_atoms = ", ".join(entry["Donor Atoms"])
                acceptor_atoms = ", ".join(entry["Acceptor Atoms"])
                row = (
                    f"{entry['Donor Index']:<12} {entry['Donor Type']:<10} {entry['Donor Orb No']:<12} "
                    f"{donor_atoms:<25} {entry['Acceptor Index']:<15} {entry['Acceptor Type']:<12} "
                    f"{entry['Acceptor Orb No']:<15} {acceptor_atoms:<25} {entry['E(2)']:>8.2f} "
                    f"{entry['E Diff']:>8.2f} {entry['Fock Elem']:>10.2f}"
                )
                print(row)
                counter += 1
        print("=" * len(header))
        print(f"Total number of LP to BD* interactions: {counter}")

##########################################
# For visualisation of NBO Second Order Perturbation Theory Analysis i am thinking of using xyz file for atom numbers with their positions and using the 
# data to draw cylinder connections between the donor and acceptor atoms with thickness and colour depending on the E(2) value
## i.e. the larger the E(2) value the thicker the cylinder and the more red it is
# i.e. the smaller the E(2) value the thinner the cylinder and the more blue it is
##########################################
    def visualise_nbo_data(self, xyz_file, view=None, display=True, donor=None, acceptor=None, donor_type="LP", acceptor_type="BD*", E2_below=None, E2_above=None, label=True, print_latex=False):
        print("*" * 150)
        print("     Note this defaults to LP to BD* interactions, if you want to see other interactions please specify the donor and acceptor types.")
        print("*" * 150)
        print("")
        connection_indexes = []
        vmin = 0
        vmax = 1
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        coordinates = []
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                coordinates.append((x, y, z))

        # Create a view for visualization
        if view is None:
            view = py3Dmol.view(width=1000, height=800)
            view.addModel(open(xyz_file, 'r').read(), 'xyz')
            view.setStyle({'stick': {'radius': 0.03}})
            view.setBackgroundColor('white')

        # Add spheres for atoms
        # for i, coord in enumerate(coordinates):
        #     view.addSphere({'center': {'x': coord[0], 'y': coord[1], 'z': coord[2]}, 'radius': 0.2, 'color': 'gray'})
        
        if donor is not None:
            if isinstance(donor, str) and len(donor) == 1:
                donor = [donor]
            if isinstance(donor, str) and len(donor) in [1, 2]:
                donor = list(donor)

        if acceptor is not None:
            if isinstance(acceptor, str) and len(acceptor) == 1:
                acceptor = [acceptor]
            if isinstance(acceptor, str) and len(acceptor) in [1, 2]:
                acceptor = list(acceptor)

        for entry in self.nbo_data:
            if entry["Donor Type"] == donor_type and entry["Acceptor Type"] == acceptor_type:
                if donor is not None:
                    if len(donor) == 1:
                        if not any(donor[0] in atom for atom in entry["Donor Atoms"]):
                            continue
                    if len(donor) == 2:
                        if len(entry["Donor Atoms"]) != 2:
                            continue
                        elif not (donor[0] == ''.join(filter(str.isalpha, entry["Donor Atoms"][0])) and donor[1] == ''.join(filter(str.isalpha, entry["Donor Atoms"][1]))):
                            continue
                if acceptor is not None:
                    if len(acceptor) == 1:
                        if not any(acceptor[0] in atom for atom in entry["Acceptor Atoms"]):
                            continue
                    if len(acceptor) == 2:
                        if len(entry["Acceptor Atoms"]) != 2:
                            continue
                        elif not (acceptor[0] == ''.join(filter(str.isalpha, entry["Acceptor Atoms"][0])) and acceptor[1] == ''.join(filter(str.isalpha, entry["Acceptor Atoms"][1]))):
                            continue
                if E2_below is not None and entry["E(2)"] > E2_below:
                    continue
                if E2_above is not None and entry["E(2)"] < E2_above:
                    continue

                donor_atoms = entry["Donor Atoms"]
                acceptor_atoms = entry["Acceptor Atoms"]
                # if donor or acceptor:
                #     if len(donor_atoms) != len(donor) or len(acceptor_atoms) != len(acceptor):
                #         continue

                donor_index = int(re.findall(r'\d+', donor_atoms[-1])[-1]) - 1
                acceptor_index = int(re.findall(r'\d+', acceptor_atoms[-1])[-1]) - 1
                e2_value = entry["E(2)"]
                norm = colors.Normalize(vmin=E2_above if E2_above else 0, vmax=E2_below if E2_below else 1)  
                cmap = plt.colormaps.get_cmap('rainbow')
                rgb = cmap(norm(e2_value))[:3]  # Extract the RGB components
                color = '#%02x%02x%02x' % (int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
                radius = 0.05

                # show e2 value as label at midpoint of the cylinder
                mid_x = (coordinates[donor_index][0] + coordinates[acceptor_index][0]) / 2
                mid_y = (coordinates[donor_index][1] + coordinates[acceptor_index][1]) / 2
                mid_z = (coordinates[donor_index][2] + coordinates[acceptor_index][2]) / 2
                if label:
                    view.addLabel(f"E(2): {e2_value:.2f}", {
                        'position': {'x': mid_x, 'y': mid_y, 'z': mid_z},
                        'backgroundColor': color,
                        'backgroundOpacity': 0.3,
                        'fontSize': 10,
                        'fontColor': 'black',
                        'fontWeight': 'bold'
                    })
                view.addCylinder({
                    'start': {'x': coordinates[donor_index][0], 'y': coordinates[donor_index][1], 'z': coordinates[donor_index][2]},
                    'end': {'x': coordinates[acceptor_index][0], 'y': coordinates[acceptor_index][1], 'z': coordinates[acceptor_index][2]},
                    'color': color,
                    'radius': radius,
                    'opacity': 0.8
                })
                connection_indexes.append((donor_index, acceptor_index))
                if entry["E(2)"] > vmax:
                    vmax = entry["E(2)"]
                if entry["E(2)"] < vmin:
                    vmin = entry["E(2)"]

        # print only the visualised data in a table
        print(f"{'Donor Index':<12} {'Donor Type':<10} {'Donor Orb No':<12} "
            f"{'Donor Atoms':<25} {'Acceptor Index':<15} {'Acceptor Type':<12} "
            f"{'Acceptor Orb No':<15} {'Acceptor Atoms':<25} {'E(2)':>8} {'E Diff':>8} {'Fock Elem':>10}")
        print("=" * 160)
        counter = 0
        for entry in self.nbo_data:
            if entry["Donor Type"] == donor_type and entry["Acceptor Type"] == acceptor_type:
                if donor is not None:
                    if len(donor) == 1:
                        if not any(donor[0] in atom for atom in entry["Donor Atoms"]):
                            continue
                    elif len(donor) == 2:
                        if len(entry["Donor Atoms"]) != 2:
                            continue
                        elif not (donor[0] == ''.join(filter(str.isalpha, entry["Donor Atoms"][0])) and donor[1] == ''.join(filter(str.isalpha, entry["Donor Atoms"][1]))):
                            continue
                if acceptor is not None:
                    if len(acceptor) == 1:
                        if not any(acceptor[0] in atom for atom in entry["Acceptor Atoms"]):
                            continue
                    elif len(acceptor) == 2:
                        if len(entry["Acceptor Atoms"]) != 2:
                            continue
                        elif not (acceptor[0] == ''.join(filter(str.isalpha, entry["Acceptor Atoms"][0])) and acceptor[1] == ''.join(filter(str.isalpha, entry["Acceptor Atoms"][1]))):
                            continue
                if E2_below is not None and entry["E(2)"] > E2_below:
                    continue
                if E2_above is not None and entry["E(2)"] < E2_above:
                    continue
                donor_atoms = ", ".join(entry["Donor Atoms"])
                acceptor_atoms = ", ".join(entry["Acceptor Atoms"])
                row = (
                    f"{entry['Donor Index']:<12} {entry['Donor Type']:<10} {entry['Donor Orb No']:<12} "
                    f"{donor_atoms:<25} {entry['Acceptor Index']:<15} {entry['Acceptor Type']:<12} "
                    f"{entry['Acceptor Orb No']:<15} {acceptor_atoms:<25} {entry['E(2)']:>8.2f} "
                    f"{entry['E Diff']:>8.2f} {entry['Fock Elem']:>10.2f}"
                )
                print(row)
                counter += 1
        print("=" * 160)
        print(f"Total number of NBO interactions of interest: {counter}")

        if print_latex:
            """prints same as above but in latex table format ready for copy and paste without the donor/acceptor index and orbital numbers"""
            print("\\begin{table}[H]")
            print("\\centering")
            print("\\begin{tabular}{|c|c|c|c|c|c|c|}")
            print("\\hline")
            print(f"{'Donor Type':<10} & {'Donor Atoms':<25} & {'Acceptor Type':<12} & "
                f"{'Acceptor Atoms':<25} & {'E(2)':>8} & {'E Diff':>8} & {'Fock Elem':>10} \\\\")
            print("\\hline")
            for entry in self.nbo_data:
                if entry["Donor Type"] == donor_type and entry["Acceptor Type"] == acceptor_type:
                    if donor is not None:
                        if len(donor) == 1:
                            if not any(donor[0] in atom for atom in entry["Donor Atoms"]):
                                continue
                        elif len(donor) == 2:
                            if len(entry["Donor Atoms"]) != 2:
                                continue
                            elif not (donor[0] == ''.join(filter(str.isalpha, entry["Donor Atoms"][0])) and
                                    donor[1] == ''.join(filter(str.isalpha, entry["Donor Atoms"][1]))):
                                continue
                    if acceptor is not None:
                        if len(acceptor) == 1:
                            if not any(acceptor[0] in atom for atom in entry["Acceptor Atoms"]):
                                continue
                        elif len(acceptor) == 2:
                            if len(entry["Acceptor Atoms"]) != 2:
                                continue
                            elif not (acceptor[0] == ''.join(filter(str.isalpha, entry["Acceptor Atoms"][0])) and
                                    acceptor[1] == ''.join(filter(str.isalpha, entry["Acceptor Atoms"][1]))):
                                continue
                    if E2_below is not None and entry["E(2)"] > E2_below:
                        continue
                    if E2_above is not None and entry["E(2)"] < E2_above:
                        continue
                    donor_atoms = ", ".join(entry["Donor Atoms"])
                    acceptor_atoms = ", ".join(entry["Acceptor Atoms"])
                    row = (f"{entry['Donor Type']:<10} & {donor_atoms:<25} & {entry['Acceptor Type']:<12} "
                        f"& {acceptor_atoms:<25} & {entry['E(2)']:>8.2f} & {entry['E Diff']:>8.2f} & {entry['Fock Elem']:>10.2f} \\\\")
                    print(row)
            print("\\hline")
            print("\\end{tabular}")
            print("\\caption{NBO Second Order Perturbation Theory Analysis}")
            print("\\label{tab:nbo_sop}")
            print("\\end{table}")
        # Show the view
        if display:
            view.zoomTo()
            view.show()
        
        # -- color bar to show the spread of E(2) values --
        fig, ax = plt.subplots(figsize=(6, 1))
        fig.subplots_adjust(bottom=0.5)
        cmap = cm.rainbow
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='horizontal')
        cb.set_label('E(2) kcal/mol range')
        plt.show()

        connection_indexes = set(connection_indexes)
        return connection_indexes if not display else None