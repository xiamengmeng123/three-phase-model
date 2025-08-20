#先激活环境 conda activate pmx
#先转换成pdb文件   gmx_mpi editconf -f complex-equil.gro -o complex-equil.pdb
#修改complex-equil.pdb中的HMX和LIG残基名
import mdtraj as md
import pmx,os


pdb_file = 'complex-equil.pdb'


# 用 pmx 读取结构文件
complex = pmx.Model(pdb_file,rename_atoms=True)

for atom in complex.atoms[:5]:
    # print(dir(atom))
    print(f"Atom id: {atom.id}, name: {atom.name}, resname: {atom.resname}, element: {atom.symbol}, mass: {atom.m}")

# 创建约束对象，ligname要对应结构中的配体残基名，如 'LIG'
rest = pmx.alchemy.AbsRestraints(complex=complex, ligname='LIG',seed=42)

# 自动选点并写 restraint 文件
rest.write_summary('restraints.info')
ii = rest.make_ii()
# # ✅ 写入 GROMACS 支持的跨分子 restraint 形式  采用6 1 9的方式
# with open('ready_top/boresch_restraints.itp', 'w') as f:
#     f.write('[ intermolecular_interactions ]\n\n')

#     # bonds
#     if 'bonds' in ii:
#         f.write('[ bonds ]\n')
#         f.write('; ai  aj  funct  length_A  force_A  length_B  force_B\n')
#         for b in ii['bonds']:
#             length_A = b[3][0]
#             force_A = b[3][3]
#             length_B = b[3][2]
#             force_B = b[3][3]
#             f.write(f"{b[0]} {b[1]} 6 {length_A:.5f} {force_A:.2f} {length_B:.5f} {force_B:.2f}\n")
#         f.write('\n')

#     # angles
#     if 'angles' in ii:
#         f.write('[ angles ]\n')
#         f.write('; ai  aj  ak  funct  angle_A  force_A  angle_B  force_B\n')
#         for a in ii['angles']:
#             angle_A = a[4][0]
#             force_A = a[4][3]
#             angle_B = a[4][2]
#             force_B = a[4][3]
#             f.write(f"{a[0]} {a[1]} {a[2]} 1 {angle_A:.3f} {force_A:.2f} {angle_B:.3f} {force_B:.2f}\n")
#         f.write('\n')

#     # dihedrals
#     if 'dihedrals' in ii:
#         f.write('[ dihedrals ]\n')
#         f.write('; ai  aj  ak  al  funct  phi_A  k_A  mult_A  phi_B  k_B  mult_B\n')
#         for d in ii['dihedrals']:
#             phi_A = d[5][0]
#             k_A = d[5][3]
#             phi_B = d[5][2]  # 可以设置成 phi_A，也可保持对称性写 phi_B
#             k_B = d[5][3]
#             multiplicity = 1
#             f.write(f"{d[0]} {d[1]} {d[2]} {d[3]} 9 {phi_A:.3f} {k_A:.2f} {multiplicity} {phi_B:.3f} {k_B:.2f} {multiplicity}\n")

#  ✅ 写入 GROMACS 支持的跨分子 restraint 形式  采用6 1 2的方式


# ✅ 写入 GROMACS 支持的跨分子 restraint 形式  采用6 1 2的方式
with open('ready_top/boresch_restraints.itp', 'w') as f:
    f.write('[ intermolecular_interactions ]\n')
    
    # 状态A（λ=0）的固定值
    STATE_A_FORCE = 0.0  # λ=0 时力常数必须为0
    STATE_A_REF = 0.0    # λ=0 时参考值设为0（实际不起作用）

    # bonds
    if 'bonds' in ii:
        f.write('[ bonds ]\n')
        f.write('; ai   aj     type  bA     kA   bB     kB\n')
        for b in ii['bonds']:
            # 状态A：参考值设为0，力常数设为0
            # 状态B：使用自动计算的参考值和力常数
            length_B = b[3][2]  # 状态B的参考键长
            force_B = b[3][3]   # 状态B的力常数
            f.write(f"{b[0]} {b[1]} 6 {STATE_A_REF:.5f} {STATE_A_FORCE:.2f} {length_B:.5f} {force_B:.2f}\n")
        f.write('\n')

    # angles
    if 'angles' in ii:
        f.write('[ angles ]\n')
        f.write('; ai   aj   ak     type  thA     kA   thB     kB\n')
        for a in ii['angles']:
            # 状态A：参考值设为0，力常数设为0
            # 状态B：使用自动计算的参考值和力常数
            angle_B = a[4][2]  # 状态B的参考角度
            force_B = a[4][3]  # 状态B的力常数
            f.write(f"{a[0]} {a[1]} {a[2]} 1 {STATE_A_REF:.3f} {STATE_A_FORCE:.2f} {angle_B:.3f} {force_B:.2f}\n")
        f.write('\n')

    # dihedrals
    if 'dihedrals' in ii:
        f.write('[ dihedrals ]\n')
        f.write('; ai   aj   ak   al      type  phiA     kA   phiB     kB\n')
        for d in ii['dihedrals']:
            # 状态A：参考值设为0，力常数设为0
            # 状态B：使用自动计算的参考值和力常数
            phi_B = d[5][2]  # 状态B的参考二面角
            k_B = d[5][3]    # 状态B的力常数
            f.write(f"{d[0]} {d[1]} {d[2]} {d[3]} 2 {STATE_A_REF:.3f} {STATE_A_FORCE:.2f} {phi_B:.3f} {k_B:.2f}\n")

# -------- 可视化辅助：打印选中的原子对对应的原子信息 --------
atom_id_to_info = {atom.id: (atom.name, atom.resname, atom.chain_id) for atom in complex.atoms}

print("\n--- Selected atoms in bonds ---")
for bond in ii.get('bonds', []):
    a1, a2 = bond[0], bond[1]
    print(f"Bond between atom {a1} {atom_id_to_info.get(a1)} and atom {a2} {atom_id_to_info.get(a2)}")

print("\n--- Selected atoms in angles ---")
for angle in ii.get('angles', []):
    a1, a2, a3 = angle[0], angle[1], angle[2]
    print(f"Angle between atoms {a1} {atom_id_to_info.get(a1)}, {a2} {atom_id_to_info.get(a2)}, {a3} {atom_id_to_info.get(a3)}")

print("\n--- Selected atoms in dihedrals ---")
for dih in ii.get('dihedrals', []):
    a1, a2, a3, a4 = dih[0], dih[1], dih[2], dih[3]
    print(f"Dihedral between atoms {a1} {atom_id_to_info.get(a1)}, {a2} {atom_id_to_info.get(a2)}, {a3} {atom_id_to_info.get(a3)}, {a4} {atom_id_to_info.get(a4)}")

# -------- 生成PDB子集文件，仅包含选中原子 --------
selected_atoms = set()
for term_type in ['bonds', 'angles', 'dihedrals']:
    for entry in ii.get(term_type, []):
        selected_atoms.update(entry[:{'bonds':2,'angles':3,'dihedrals':4}[term_type]])

with open(pdb_file) as fin, open('selected_atoms.pdb', 'w') as fout:
    for line in fin:
        if line.startswith(('ATOM', 'HETATM')):
            atom_id = int(line[6:11].strip())
            if atom_id in selected_atoms:
                fout.write(line)

print(f"\nExtracted selected atoms to 'selected_atoms.pdb', you can visualize this file in PyMOL/VMD/etc.")

