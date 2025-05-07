### This script combines QTAIM and NBO data visualization using py3Dmol.
### It requires both QTAIM and NBO scripts to be available:
#
# from qtaim import QTAIM
# from nbo import NBO_SOP
# import py3Dmol

def combine_qtaim_and_nbo_data(qtaim_file, nbo_file, xyz_file, highlight_connections=True, donor=None, acceptor=None, donor_type="LP", acceptor_type="BD*", A=None, B=None, E2_below=None, label_nbo=True, label_qtaim=True, legend=True, print_latex=False, stickThickness=0.03, bond_lengths=False):
    qtaim = QTAIM(qtaim_file)
    nbo = NBO_SOP(nbo_file)
    view = py3Dmol.view(width=1500, height=1500)
    xyz_data = open(xyz_file, 'r').read()
    view.addModel(xyz_data, 'xyz')
    view.setStyle({'stick': {'radius': stickThickness}})
    view.setViewStyle({'style': 'perspective'})

    nbo_indexes = nbo.visualise_nbo_data(xyz_file, view=view, display=False, donor = donor, acceptor = acceptor, donor_type=donor_type, acceptor_type=acceptor_type, E2_below=E2_below, label=label_nbo, print_latex=print_latex)
    qtaim_indexes = qtaim.visualise(xyz_file, view=view, display=False,  show_pos_lap=label_qtaim, connect_atoms_A_B=True, A=A if A else donor, B=B if B else acceptor, print_parameters=True, covalent=False, legend=legend, print_latex=print_latex, show_bond_lengths=bond_lengths)

    if highlight_connections:
        specific_atoms = []
        combined_indexes = qtaim_indexes | nbo_indexes
        for tpl in combined_indexes:
            for index in tpl:
                if index not in specific_atoms:
                    specific_atoms.append(index)
        for i, line in enumerate(xyz_data.splitlines()[2:]):
            if i in specific_atoms:
                parts = line.split()
                if len(parts) >= 4:
                    x, y, z = map(float, parts[1:4])
                    color = 'blue' if parts[0] == 'N' else 'white' if parts[0] == 'H' else 'red' if parts[0] == 'O' else 'gray'
                    view.addSphere({'center': {'x': x, 'y': y, 'z': z}, 'radius': 0.2, 'color': color, 'opacity': 1})
    view.zoomTo()
    view.show()