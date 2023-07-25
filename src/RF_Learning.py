import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from deepchem.utils.vina_utils import prepare_inputs
import deepchem as dc
from rdkit import Chem
import scipy

from deepchem.utils.evaluate import Evaluator
from sklearn.ensemble import RandomForestRegressor
structure_df = pd.read_csv("../data/structures.csv")
structures = np.array(structure_df['smiles'])

affinity_df = pd.read_csv("../data/test.csv")
affinities = np.array(affinity_df[' G'])
affinities_normalized = affinities / np.linalg.norm(affinities)