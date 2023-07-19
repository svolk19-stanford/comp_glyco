def prepare_ligands():
    for i in range(94):
        cmd.load(f'pymol_fragments/out{i + 1}.sdf')
        cmd.label(f'out{i+1}', '"%s"%(index)')
        cmd.save(f'/Applications/PyMOL.app/Contents/share/pymol/data/chempy/fragments/lig_{i}.pkl', f'out{i+1}')

def assemble_ligands():
    h_attachment_pos = [18, 12, 16, 15, 18, 17, 18, 18, 16, 20, 18, 19, 18, 15, 16, 16, 19, 16, 19, 15, 19, 20, 15, 19, 18, 22, 20, 16, 19, 20, 22, 21, 14, 17, 34, 15, 18, 22, 17, 17, 16, 18, 18, 15, 17, 18, 16, 28, 15, 14, 14, 19, 15, 13, 18, 16, 15, 16, 22, 17, 19, 15, 16, 17, 24, 17, 18, 16, 18, 19, 19, 19, 14, 19, 20, 20, 20, 25, 19, 23, 21, 27, 21, 21, 27, 21, 21, 21, 19, 26, 25, 23, 24, 16]
    for i in range(len(h_attachment_pos)):
        print(i)
        cmd.load('scaffold_2g5r.pdb')
        cmd.select('pk1', '/scaffold_2g5r/C/A/NXD`145/H01')
        editor.attach_fragment('pk1',f'lig_{i}', h_attachment_pos[i] - 1,0)
        cmd.save(f'pymol_assembled_ligands/lig_assembled_{i}.pdb', 'scaffold_2g5r')
        cmd.delete('scaffold_2g5r')

cmd.extend("prepare_ligands", prepare_ligands)
cmd.extend("assemble_ligands", assemble_ligands)





