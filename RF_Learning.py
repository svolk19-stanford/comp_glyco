import pandas as pd
import numpy as np
from sklearn.preprocessing import normalize
from deepchem.utils.vina_utils import prepare_inputs
from pdbfixer import PDBFixer

structure_df = pd.read_csv("../data/structures.csv")
structures = np.array(structure_df['smiles'])
affinity_df = pd.read_csv("../data/test.csv")
affinities = np.array(affinity_df[' G'])
affinities_normalized = affinities / np.linalg.norm(affinities)

pdbids = ['1O7S' for i in range(len(structures))]
for (pdbid, ligand) in zip(pdbids, structures):
    #   fixer = PDBFixer(url='https://files.rcsb.org/download/%s.pdb' % (pdbid))
    #   PDBFile.writeFile(fixer.topology, fixer.positions, open('%s.pdb' % (pdbid), 'w'))
    p, m = prepare_inputs('%s.pdb' % (pdbid),
                          ligand,
                          replace_nonstandard_residues=False,
                          remove_heterogens=False,
                          remove_water=False,
                          add_hydrogens=False)

    if p and m:  # protein and molecule are readable by RDKit
        Chem.rdmolfiles.MolToPDBFile(p, '%s.pdb' % (pdbid))
        Chem.rdmolfiles.MolToPDBFile(m, 'ligand_%s.pdb' % (pdbid))
