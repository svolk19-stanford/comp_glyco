import os 
os.system("vina --receptor 2hrl.pdbqt --ligand 1.pdbqt \
          --config ../vina_test_data/config.txt \
          --exhaustiveness=32 \
          --local_only\
          --out ../vina_test_data/1_out.pdbqt")