import vina
import pandas as pd 
import meeko
from tqdm import tqdm

num_poses = 20
energies_optimized = []
energies = []
for i in tqdm(range(94)):
    v = vina.Vina(sf_name='vina', verbosity=0)
    v.set_receptor('../data/2g5r.pdbqt')
    v.set_ligand_from_file(f'../data/pymol_ligands/lig_assembled_{i}.pdbqt')
    # for 2hrl
    # v.compute_vina_maps(center=[24.858, 30.832, 14.323], box_size=[50, 50, 40])
    # for 2g5r
    v.compute_vina_maps(center=[0.663, 9.447, 23.69], box_size=[44, 40, 70])

    # Minimized locally the current pose
    energy_minimized = v.optimize()
    energies_optimized.append(list(energy_minimized)[0])
    output_energies_o = pd.DataFrame({'energy': energies_optimized})
    output_energies_o.to_csv("../out/pymol/output_energies_optimized_2g5r_1.csv")
    v.write_pose(f'../out/pymol/optimized_2g5r_1/ligand_{i}.pdbqt', overwrite=True)

    # Dock ligand and calculate energies
    v.dock(exhaustiveness=32, n_poses=20)
    output_pdbqt = v.poses(n_poses=20)
    energy = v.energies(num_poses)
    energy = [i[0] for i in energy]
    energies.append(energy)
    output_energies = pd.DataFrame(energies, columns=[str(i) for i in range(num_poses)])
    output_energies.to_csv("../out/pymol/output_energies_docked_2g5r_1.csv")

    # Save docked poses
    pmol = meeko.PDBQTMolecule(output_pdbqt)
    for j, pose in enumerate(pmol):
        pmol.write_pdbqt_file(f'../out/pymol/optimized_and_docked_2g5r_1/ligand_{i}_pose_{j}.pdbqt', overwrite=True)

    






