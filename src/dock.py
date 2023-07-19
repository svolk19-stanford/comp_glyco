import scipy.stats as stats
from sklearn.metrics import r2_score 
import vina
import pandas as pd 
import meeko
from rdkit import Chem
from tqdm import tqdm

affinities_df = pd.read_csv('../data/test.csv')
affinities = list(affinities_df['G'])

num_poses = 20
ligands = pd.read_csv("../data/substituted_sia.csv")
scores = []
energies = []
# for i in range(82, len(ligands['smiles'])):
#     lig = ligands['smiles'][i]
for i in tqdm(range(65, len(ligands['smiles']))):
    lig = ligands['smiles'][i]
    mol = Chem.MolFromSmiles(lig)
    protonated_lig = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(protonated_lig)
    meeko_prep = meeko.MoleculePreparation()
    meeko_prep.prepare(protonated_lig)
    lig_pdbqt = meeko_prep.write_pdbqt_string()
    v = vina.Vina(sf_name='vina', verbosity=0)
    v.set_receptor('../data/2hrl_rigid_multiple.pdbqt', '../data/2hrl_flex_multiple.pdbqt')
    v.set_ligand_from_string(lig_pdbqt)
    # from gnina
    v.compute_vina_maps(center=[24.839, 30.851, 14.8775], box_size=[22.38, 26.076, 21.429])

    # bounding box for 2df3
    # v.compute_vina_maps(center=[18.192, 34.378, 23.055], box_size=[42, 42, 48])
    # bounding box for 2hrl: semi-flex
    # v.compute_vina_maps(center=[26.817, 33.939, 17.122], box_size=[28, 28, 28])
    # bounding box for 2hrl: fully flex
    # v.compute_vina_maps(center=[28.409, 32.3, 14.201], box_size=[26, 38, 24])
    # # bounding box for 2hrl: fragment region only
    # v.compute_vina_maps(center=[32.395, 34.447, 10.547], box_size=[24, 28, 36])
    # v.optimize()
    v.dock(exhaustiveness=32, n_poses=20)
    output_pdbqt = v.poses(n_poses=20)
    energy = v.energies(num_poses)
    # print(energy)
    # assert(False)
    energy = [i[0] for i in energy]
    energies.append(energy)
    output_energies = pd.DataFrame(energies, columns=[str(i) for i in range(num_poses)])
    output_energies.to_csv("../out/2hrl_poses_subst_sia/output_energies.csv")

    pmol = meeko.PDBQTMolecule(output_pdbqt)
    for j, pose in enumerate(pmol):
        pmol.write_pdbqt_file(f'../out/2hrl_poses_subst_sia/ligand_{i}_pose_{j}_2HRL.pdbqt', overwrite=True)

    assert(False)
    # y_hat = [j[0] for j in energies]
    # y = affinities[:(i + 1)]
    # slope, intercept, r_value, p_value, std_err = stats.linregress(y,y_hat)
    # print(f"r value: {r_value}")
    # print(f"r^2: {r_value**2}")
    # print(f"p val: {p_value}")






