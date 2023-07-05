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
# for i in range(82, len(ligands['smiles'])):
#     lig = ligands['smiles'][i]
for i, lig in enumerate(ligands['smiles']): 
    mol = Chem.MolFromSmiles(lig)
    protonated_lig = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(protonated_lig)
    meeko_prep = meeko.MoleculePreparation()
    meeko_prep.prepare(protonated_lig)
    lig_pdbqt = meeko_prep.write_pdbqt_file(f"../out/ligand_pdbqt_sia/ligand_{i + 1}.pdbqt")
    print(f"prepared ligand {i}")
    # v = vina.Vina(sf_name='vina', verbosity=0)
    # v.set_receptor('../data/2hrl_rigid_fragment.pdbqt', '../data/2hrl_flex_fragment.pdbqt')
    # v.set_ligand_from_string(lig_pdbqt)
    # bounding box for 2df3
    # v.compute_vina_maps(center=[18.192, 34.378, 23.055], box_size=[42, 42, 48])
    # bounding box for 2hrl: semi-flex
    # v.compute_vina_maps(center=[26.817, 33.939, 17.122], box_size=[28, 28, 28])
    # bounding box for 2hrl: fully flex
    # v.compute_vina_maps(center=[28.409, 32.3, 14.201], box_size=[26, 38, 24])
    # # bounding box for 2hrl: fragment region only
    # v.compute_vina_maps(center=[32.395, 34.447, 10.547], box_size=[24, 28, 36])
    # v.optimize()
    # v.dock(exhaustiveness=32, n_poses=20)
    # output_pdbqt = v.poses(n_poses=20)
    # energy = v.energies(num_poses)
    # print(energy)
    # assert(False)
    # energy = [i[0] for i in energy]
    # energies.append(energy)
    # output_energies = pd.DataFrame(energies, columns=[str(i) for i in range(num_poses)])
    # output_energies.to_csv("../out/output_energies_2HRL_fragment_32.csv")

    # pmol = meeko.PDBQTMolecule(output_pdbqt)
    # score = []
    # for j, pose in enumerate(pmol):
    #     score.append(pmol.score)
    #     pmol.write_pdbqt_file(f'../out/2hrl_poses_fragment_32/ligand_{i}_pose_{j}_2HRL.pdbqt', overwrite=True)
    # scores.append(score)
    # output_scores = pd.DataFrame(scores, columns=[str(j) for j in range(num_poses)])
    # output_scores.to_csv("../out/output_scores_2HRL_fragment_32.csv")

    # y_hat = [j[0] for j in energies]
    # y = affinities[:(i + 1)]
    # slope, intercept, r_value, p_value, std_err = stats.linregress(y,y_hat)
    # print(f"r value: {r_value}")
    # print(f"r^2: {r_value**2}")
    # print(f"p val: {p_value}")






