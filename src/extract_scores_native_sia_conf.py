import pandas as pd
import subprocess as sp

y_hat = pd.read_csv('../out/pymol/output_energies_docked_2g5r_2.csv')
y = pd.read_csv('../data/test.csv')['G']
y_native_conf = []
yhat_native_conf = []

available_ligands = sp.check_output("cd ../out/pymol/optimized_and_docked_2g5r_2/native_sia_conf_subset && ls", shell=True)
available_ligands = available_ligands.splitlines()
print(available_ligands)

for lig in available_ligands:
    lig_split = str(lig).split("_")
    lig_num = lig_split[1]
    pose_num = lig_split[3].split('.')[0]
    y_native_conf.append(y[int(lig_num)])
    yhat_native_conf.append(y_hat.loc[int(lig_num)][pose_num])

