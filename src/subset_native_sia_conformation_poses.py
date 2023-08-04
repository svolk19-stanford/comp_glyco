import subprocess
import os

num_poses = 20
num_ligands = 94
data_path = '../out/pymol/optimized_and_docked_2g5r_2/'

for i in range(16, num_ligands):
    print(i)
    for j in range(num_poses):
        try:
            out = subprocess.check_output("obrms -f " + data_path + "scaff.pdb " + data_path + f"ligand_{i}_pose_{j}.pdbqt", shell=True)
            rmsd = float(out.split()[-1])
            if rmsd <= 2:
                os.system("cp " + data_path + f"ligand_{i}_pose_{j}.pdbqt " + data_path + f"/native_sia_conf_subset/ligand_{i}_pose_{j}.pdbqt")
        except Exception as error:
            print(error)
            continue



        