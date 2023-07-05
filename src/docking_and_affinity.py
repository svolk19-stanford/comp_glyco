import deepchem
import os
import numpy as np
import pandas as pd

import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem
import deepchem as dc

from deepchem.utils import download_url, load_from_disk


from openmm.app import PDBFile
from pdbfixer import PDBFixer

from deepchem.utils.vina_utils import prepare_inputs

pdbid = "1o7s"
ligand = "OC[C@H]1O[C@@H](OCCN=[N+]=[N-])[C@H](O)C(O)[C@@H]1O[C@@H]2O[C@H](CO[C@]3(C(O)=O)C[C@H](O)[C@@H](NC(C)=O)[C@H](C(O)C(O)CNC(OCC4=CN(c5c(C(=O)N)cccc5)N=N4)=O)O3)[C@H](O)[C@H](O)[C@H]2O"

fixer = PDBFixer(pdbid=pdbid)
PDBFile.writeFile(fixer.topology, fixer.positions, open('%s.pdb' % (pdbid), 'w'))

p, m = None, None
# fix protein, optimize ligand geometry, and sanitize molecules
try:
    p, m = prepare_inputs('%s.pdb' % (pdbid), ligand)
except:
    print('%s failed PDB fixing' % (pdbid)) 

if p and m:  # protein and molecule are readable by RDKit
    print(pdbid, p.GetNumAtoms())
    Chem.rdmolfiles.MolToPDBFile(p, '%s.pdb' % (pdbid))
    Chem.rdmolfiles.MolToPDBFile(m, 'ligand_%s.pdb' % (pdbid))

finder = dc.dock.binding_pocket.ConvexHullPocketFinder()
pockets = finder.find_pockets('1o7s.pdb')
vpg = dc.dock.pose_generation.VinaPoseGenerator()
complexes, scores = vpg.generate_poses(molecular_complex=('1o7s.pdb', 'ligand_1o7s.pdb'),  # protein-ligand files for docking,
                                       out_dir='vina_test',
                                       generate_scores=True
                                      )

print(scores)
