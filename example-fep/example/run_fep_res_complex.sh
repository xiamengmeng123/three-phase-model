#!/bin/bash 

##########
#by mmx 20240620
#yhbatch -N 2 -p cp2 run_fep_res_complex.sh 1 cp2 
##########
num=$1
partition=$2
ff=1four #提供带数字的文件夹名

FEP_dir="/mmx/gromacs"
base_dir="${FEP_dir}/NPBA/${num}/${ff}"   
ready_top_dir="${base_dir}/ready_top"
prep_dir="${base_dir}/prep"
LIG="LIG_${num}"
HMX="HMX_${num}"

command1="yhrun -N 1 -p $partition -J $num"
command2="yhrun -N 2 -p $partition -n  72 -J $num"

traj_file="${prep_dir}/md.xtc"
if [ -f "$traj_file" ]; then
    echo "预平衡轨迹文件已存在"
else
    mkdir -p "${prep_dir}"
    cd $prep_dir
    echo "$ready_top_dir/complex.gro"
    cp -rf "${ready_top_dir}/complex.gro" $prep_dir
    cp -rf "${ready_top_dir}/forcefield.itp" $prep_dir
    cp -rf "${ready_top_dir}/topol.top" $prep_dir
    cp -rf "${ready_top_dir}/$LIG.gro" $prep_dir
    cp -rf "${ready_top_dir}/$HMX.gro" $prep_dir
    cp -rf "${ready_top_dir}/$LIG.itp" $prep_dir
    cp -rf "${ready_top_dir}/$HMX.itp" $prep_dir
    cp -rf ${FEP_dir}/mdp/*.mdp $prep_dir

    sed -i  "28s/.*/couple\-moltype          \= ${LIG}/"  $prep_dir/em.mdp
    sed -i "63s/.*/couple\-moltype           \= ${LIG}/"  $prep_dir/nvt.mdp
    sed -i "57s/.*/couple\-moltype           \= ${LIG}/"  $prep_dir/md.mdp

    #替换每个.gro文件中包含#include行（7，8行）的路径
    sed -i  "7s/.*/#include \"${ready_top_dir//\//\\/}\/forcefield.itp\"/" topol.top
    # sed -i  "8s/.*/#include \"${ready_top_dir//\//\\/}\/$HMX.itp\"/" topol.top
    # sed -i  "9s/.*/#include \"${ready_top_dir//\//\\/}\/$LIG.itp\"/" topol.top
    awk -v hmx_path="$ready_top_dir/$HMX.itp" \
    -v lig_path="$ready_top_dir/$LIG.itp" \
    -v hmx_posre="$prep_dir/posre_HMX.itp" \
    -v lig_posre="$prep_dir/posre_LIG.itp" \
'{
    if (NR == 8) {
        print "#include \"" hmx_path "\""
        print "; Position restraints for HMX"
        print "#ifdef POSRES"
        print "#include \""hmx_posre "\""
        print "#endif\n"
    }
    else if (NR == 9) {
        print "#include \"" lig_path "\""
        print "; Position restraints for Ligand"
        print "#ifdef POSRES"
        print "#include \"" lig_posre "\""
        print "#endif"
    }
    else {
        print
    }
}' topol.top > temp.top && mv temp.top topol.top

###########
#预平衡
###########
    echo "=====================开始预平衡==================="
    cp -ir $FEP_dir/ions.mdp $prep_dir
    ${command1} gmx_mpi editconf -f complex.gro -o newbox.gro  -bt triclinic -box 6.6414812 3.9240 8 -angles 90 90 82.91943 -noc
    # ${command1} gmx_mpi solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solvated.gro
    # ${command1} gmx_mpi grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
    # ${command1} gmx_mpi genion -s ions.tpr -p topol.top -o solvated_ions.gro -neutral -conc 0.15
    
    ${command1} gmx_mpi make_ndx -f $HMX.gro -o index_HMX_heavy.ndx <<EOF
! a H*
name 3 non-Hydrogen
q
EOF

    ${command1} gmx_mpi make_ndx -f $LIG.gro -o index_LIG_heavy.ndx <<EOF
! a H*
name 3 non-Hydrogen
q
EOF

    echo "non-Hydrogen" |${command1} gmx_mpi genrestr -f $HMX.gro -n index_HMX_heavy.ndx -o posre_HMX.itp
    echo "non-Hydrogen" |${command1} gmx_mpi genrestr -f $LIG.gro -n index_LIG_heavy.ndx -o posre_LIG.itp

    

    #em
    echo "=========开始能量最小化=============="
    export OMP_NUM_THREADS=1 
    ${command1} gmx_mpi grompp -f em.mdp -c newbox.gro  -p topol.top  -o em.tpr -maxwarn 1 
    ${command2} gmx_mpi mdrun -v -deffnm em 
    #NVT
    echo "=========开始NVT模拟=============="   
    ${command1} gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro  -p topol.top  -o nvt.tpr -maxwarn 1   
    ${command2} gmx_mpi mdrun -v -deffnm nvt   
    # #NPT
    # echo "=========开始NPT模拟=============="   
    # ${command1} gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top  -o npt.tpr -maxwarn 1   
    # ${command2} gmx_mpi mdrun -v -deffnm npt   
    #MD
    echo "=========开始MD模拟=============="
    ${command1} gmx_mpi grompp -f md.mdp -c nvt.gro   -p topol.top  -o md.tpr -maxwarn 1   
    ${command2} gmx_mpi mdrun -v -deffnm md   
    # echo "0" |${command1} gmx_mpi trjconv -f md.xtc -s md.tpr -pbc mol -ur compact 
  
    ${command1} gmx_mpi make_ndx -f md.tpr -o index.ndx <<EOF
! a 1-8400
name 3 LIG
a 1-8400
name 4 HMX
q
EOF

    echo -e "LIG\nSystem" |${command1} gmx_mpi trjconv -f md.xtc -s md.tpr -n index.ndx -o lastframe.gro -dump -1 -pbc mol -ur compact -center 
    cp -rf lastframe.gro "${base_dir}/complex-equil.gro"
    
    rm -rf \#*

    # cp -rf md.gro "${base_dir}/complex-equil.gro"
    cp -rf topol.top $base_dir/complex-res.top    
    echo "; include restraint 文件（必须放在 molecules 段之后）" >> $base_dir/complex-res.top
    echo "#include \"${ready_top_dir}/boresch_restraints.itp\"" >> $base_dir/complex-res.top 
       
fi


################
#准备微扰文件
################
cd $base_dir
cp -rf $FEP_dir/select_boresch.py ./
${command1} gmx_mpi editconf -f complex-equil.gro -o complex-equil.pdb
awk '{
  # 处理 ATOM 行
  if (substr($0,1,4) == "ATOM") {
    # 提取原子序号 (第7-11列)
    serial = substr($0,7,5) + 0  # +0 转换为数字
    
    # 确定残基名称
    resname = (serial <= 8400) ? "HMX" : "LIG"
    
    # 替换残基名称 (第18-20列)
    printf "%s%-3s%s\n", 
           substr($0,1,17), 
           resname, 
           substr($0,21)
  } 
  # 处理其他行 (TER/ENDMDL等)
  else {
    print
  }
}' complex-equil.pdb > modified.pdb   && mv modified.pdb complex-equil.pdb
${command1} python select_boresch.py

mkdir -p  MDP 
cp -rf  $FEP_dir/MDP/em_steep.mdp $base_dir/MDP
cp -rf  $FEP_dir/MDP/nvt.mdp $base_dir/MDP
cp -rf  $FEP_dir/MDP/md.mdp $base_dir/MDP

sed -i  "43s/.*/couple\-moltype          \= ${LIG}/"  $base_dir/MDP/em_steep.mdp
# sed -i "42s/LIG/${LIG}/g" em_steep.mdp
sed -i "54s/.*/couple\-moltype           \= ${LIG}/"  $base_dir/MDP/nvt.mdp
sed -i "53s/.*/couple\-moltype           \= ${LIG}/"  $base_dir/MDP/md.mdp

cp -rf $FEP_dir/MDP/write_mdp.pl  $base_dir/MDP
cd $base_dir/MDP
num_lambda=27
perl write_mdp.pl  em_steep.mdp  $num_lambda && perl write_mdp.pl  nvt.mdp  $num_lambda   && perl write_mdp.pl  md.mdp  $num_lambda
cd $base_dir

###############
#自由能微扰
###############
# cp -rf $FEP_dir/GMX-FEP-com.sh $base_dir
# sed -i "32s/.*/loops\=${num_lambda}/" $base_dir/GMX-FEP-com.sh
# # yhbatch  -p $partition -J $num -o $base_dir/slurm.out  GMX-FEP-com.sh  $partition
# sh  GMX-FEP-com.sh  $partition

cp -rf $FEP_dir/GMX-FEP-com-lambda.sh $base_dir
#生成作业列表
rm -rf \#*
rm -rf job_list.txt
for i in `seq 0 $num_lambda`;do echo "yhbatch -N 2  -p $partition -J ${ff}_${i} -o $base_dir/slurm_%j.out GMX-FEP-com-lambda.sh $partition $i $base_dir" >> job_list.txt;done
# # 分批提交作业
# submitted=0
# while [ $submitted -lt $num_lambda ]; do
#     # 检查当前作业数量
#     current_jobs=$(yhq -u chengkun_wu |grep $partition |wc -l)
    
#     # 计算可提交数量
#     available_slots=$((15 - current_jobs))   #一个分区最多15个作业
#     if [ $available_slots -gt 0 ]; then
#         # 计算实际可提交数量（不超过剩余作业数）
#         to_submit=$((available_slots < (num_lambda - submitted) ? available_slots : (num_lambda - submitted)))
        
#         # 提交作业
#         start_index=$((submitted + 1))
#         end_index=$((submitted + to_submit))
        
#         echo "提交作业 $submitted 到 $((submitted + to_submit - 1)) (可用槽位: $available_slots, 提交数量: $to_submit)"
#         sed -n "${start_index},${end_index}p" job_list.txt | sh
        
#         submitted=$((submitted + to_submit))
#     else
#         # 等待可用槽位
#         echo "当前作业已达上限 ($current_jobs/$max_jobs), 等待 $check_interval 秒..."
#         sleep $check_interval
#     fi
# done

# echo "所有 $num_lambda 个作业已提交"

# #############
# #bar处理
# #############
# mkdir -p $base_dir/bar
# for i in `seq 0 $num_lambda`;do cp -rf Lambda_$i/Production_MD/md$i.xvg bar;done
# cd $base_dir/bar
# ${command1}  gmx_mpi bar -f md{0..41}.xvg -o -oi -oh

cd $FEP_dir
