import scipy.stats as stats
from sklearn.metrics import r2_score 
import vina
import pandas as pd 
import meeko
from rdkit import Chem
from tqdm import tqdm


for i in tqdm(range(94)):
    mol = Chem.MolFromPDBFile(f"../data/pymol_assembled_ligands/lig_assembled_{i}.pdb")
    protonated_lig = Chem.AddHs(mol)
    Chem.MolToPDBFile(protonated_lig, f"../data/pymol_ligands/lig_{i}.pdbqt")
    del protonated_lig, mol

