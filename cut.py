import numpy as np
import sys
import sh

n_atoms = 500500
n_mols = 500
header_lines = 9
lmp_file = "cnt.lmp"
data_file = "data.iso_cut"

lines = list(sh.tail(f"-{n_atoms + header_lines}", lmp_file))
header = lines[:9]
data = lines[9:]
data = np.array([[float(word) for word in line.split()] for line in data])
data = data[np.argsort(data[:, 0])]

lim = np.array([[float(word) for word in line.split()] for line in header[5:8]])
half_box = 0.5 * (lim[:, 1] - lim[:, 0])

mols = []

for i in range(1, n_mols+1):
  mol_data = data[data[:, 2] == i]
  n_atoms_mol = len(mol_data)

  mols += [[]]
  
  for j in range(n_atoms_mol - 1):
    mols[-1] += [mol_data[j]]

    x_cut = abs(mol_data[j + 1, 3] - mol_data[j, 3]) > half_box[0]
    y_cut = abs(mol_data[j + 1, 4] - mol_data[j, 4]) > half_box[1]
    if x_cut or y_cut:
      mols += [[]]
    
  mols[-1] += [mol_data[-1]]

mols = [mol for mol in mols if len(mol) > 10]

atoms = []
bonds = []
angles = []

atom_id = 1
bond_id = 1
angle_id = 1

for i, mol in enumerate(mols):
  
  atoms += [[atom_id, i+1, 1, 0, mol[0][3], mol[0][4], mol[0][5]]]
  bonds += [[bond_id, 1, atom_id, atom_id+1]]
  angles += [[angle_id, 1, atom_id, atom_id+1, atom_id+2]]
  atom_id += 1
  bond_id += 1
  angle_id += 1

  for atom in mol[1:-2]:
    atoms += [[atom_id, i+1, 2, 0, atom[3], atom[4], atom[5]]]
    bonds += [[bond_id, 1, atom_id, atom_id+1]]
    angles += [[angle_id, 1, atom_id, atom_id+1, atom_id+2]]
    atom_id += 1
    bond_id += 1
    angle_id += 1
    
  atoms += [[atom_id, i+1, 2, 0, mol[-2][3], mol[-2][4], mol[-2][5]]]
  bonds += [[bond_id, 1, atom_id, atom_id+1]]
  atom_id += 1
  bond_id += 1
  
  atoms += [[atom_id, i+1, 1, 0, mol[-1][3], mol[-1][4], mol[-1][5]]]
  atom_id += 1

with open(data_file, 'w') as f:
  f.write("LAMMPS CNT Data File\n\n")
  
  f.write(f"  {len(atoms)} atoms\n")
  f.write(f"  {len(bonds)} bonds\n")
  f.write(f"  {len(angles)} angles\n")
  f.write("  0 dihedrals\n")
  f.write("  0 impropers\n\n")
  
  f.write("  2 atom types\n")
  f.write("  1 bond types\n")
  f.write("  1 angle types\n")
  f.write(f"  {lim[0, 0]} {lim[0, 1]} xlo xhi\n")
  f.write(f"  {lim[1, 0]} {lim[1, 1]} ylo yhi\n")
  f.write(f"  {lim[2, 0]} {lim[2, 1]} zlo zhi\n\n")
  
  f.write("Masses\n\n")

  f.write("  1 975.117\n")
  f.write("  2 1950.23\n\n")

  f.write("Atoms\n\n")

  for atom in atoms:
    for entry in atom:
      f.write(f"{entry} ")
    f.write("\n")
  
  f.write("\nBonds\n\n")

  for bond in bonds:
    for entry in bond:
      f.write(f"{entry} ")
    f.write("\n")
  
  f.write("\nAngles\n\n")

  for angle in angles:
    for entry in angle:
      f.write(f"{entry} ")
    f.write("\n")
