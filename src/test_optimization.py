import scipy.stats as stats
from sklearn.metrics import r2_score 
import vina
import pandas as pd 
import meeko
from rdkit import Chem

affinities_df = pd.read_csv('../data/test.csv')
affinities = list(affinities_df['G'])

num_poses = 20
ligands = pd.read_csv("../data/substituted_sia.csv")
scores = []
energies = []
to_test = [52, 1, 32, 41, 28, 48, 31, 19, 35]
for i, lig in enumerate(ligands['smiles']): 
    if i + 1 not in to_test:
        continue
    mol = Chem.MolFromSmiles(lig)
    protonated_lig = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(protonated_lig)
    meeko_prep = meeko.MoleculePreparation()
    meeko_prep.prepare(protonated_lig)
    lig_pdbqt = meeko_prep.write_pdbqt_string()
    print(f"prepared ligand {i}")
    v = vina.Vina(sf_name='vina', verbosity=0)
    v.set_receptor('../data/2hrl_rigid_multiple.pdbqt', '../data/2hrl_flex_multiple.pdbqt')
    v.set_ligand_from_string(lig_pdbqt)
    # bounding box for 2df3
    # v.compute_vina_maps(center=[18.192, 34.378, 23.055], box_size=[42, 42, 48])
    # bounding box for 2hrl: semi-flex
    # v.compute_vina_maps(center=[26.817, 33.939, 17.122], box_size=[28, 28, 28])
    # bounding box for 2hrl: fully flex
    v.compute_vina_maps(center=[28.409, 32.3, 14.201], box_size=[26, 38, 24])
    # # bounding box for 2hrl: fragment region only
    # v.compute_vina_maps(center=[32.395, 34.447, 10.547], box_size=[24, 28, 36])
    v.dock(exhaustiveness=32, n_poses=20)
    # v.optimize()
    output_pdbqt = v.poses(n_poses=20)
    energy = v.energies(num_poses)
    print(energies)
    energy = [i[0] for i in energy].append(affinities[i])
    energies.append(energy)
    output_energies = pd.DataFrame(energies, columns=[str(i) for i in range(num_poses)].append("y"))
    output_energies.to_csv("../out/energies_2HRL_sia_optimized.csv")

    pmol = meeko.PDBQTMolecule(output_pdbqt)
    for j, pose in enumerate(pmol):
        pmol.write_pdbqt_file(f'../out/2hrl_poses_sia_optimized/ligand_{i}_pose_{j}_2HRL.pdbqt', overwrite=True)







