import scipy.stats as stats
from sklearn.metrics import r2_score 
import vina
import pandas as pd 
import meeko
from rdkit import Chem
from tqdm import tqdm


num_poses = 20
v = vina.Vina(sf_name='vina', verbosity=0)
v.set_receptor('../data/2hrl_rigid_multiple.pdbqt', '../data/2hrl_flex_multiple.pdbqt')
v.set_ligand_from_file('../data/lig_assembled.pdbqt')
v.compute_vina_maps(center=[24.839, 30.851, 14.8775], box_size=[22.38, 26.076, 21.429])

# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# # Minimized locally the current pose
# energy_minimized = v.optimize()
# print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
# v.write_pose('../out/lig_assembled_minimized.pdbqt', overwrite=True)

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=20)
v.write_poses('../out/lig_assembled_vina_out_unmin.pdbqt', n_poses=5, overwrite=True)







