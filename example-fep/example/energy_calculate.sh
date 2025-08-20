#!/bin/bash
# yhbatch -N 1 -p cp1  energy_calculate.sh  cp1

command1="yhrun -N 1 -p $1"
command2="yhrun -N 1 -p $1 -n  36"



energy_HMX(){
 
    ${command1} gmx_mpi make_ndx -f md.tpr -o HMX.ndx << EOF 
a 1-8400
name 3 HMX
q
EOF
    echo 3 |${command1} gmx_mpi trjconv -f md.xtc -s md.tpr -n HMX.ndx -o HMX.xtc  
    echo 3 |${command1} gmx_mpi convert-tpr -s md.tpr -n HMX.ndx  -o HMX.tpr  
    ${command2} gmx_mpi mdrun -s HMX.tpr  -rerun HMX.xtc -e HMX.edr 
    ${command1} gmx_mpi energy -f HMX.edr -o HMX.xvg  << EOF
LJ-(SR)
Disper.-corr.
Coulomb-(SR)
Coul.-recip.
EOF
}

energy_ligand(){
 
    ${command1} gmx_mpi make_ndx -f md.tpr -o ligand.ndx << EOF 
! a 1-8400
name 3 ligand
q
EOF
    echo 3 |${command1} gmx_mpi trjconv -f md.xtc -s md.tpr -n ligand.ndx -o ligand.xtc   
    echo 3 |${command1} gmx_mpi convert-tpr -s md.tpr -n ligand.ndx  -o ligand.tpr   
    ${command2} gmx_mpi mdrun -s ligand.tpr  -rerun ligand.xtc -e ligand.edr  
    ${command1} gmx_mpi energy -f ligand.edr -o ligand.xvg   << EOF
LJ-(SR)
Disper.-corr.
Coulomb-(SR)
Coul.-recip.
EOF
} 


    

#提取complex、HMX、ligand能量，得到三个xvg文件 
echo "=========开始计算相互作用能======="
    ${command1} gmx_mpi energy -f md.edr -o complex.xvg   << EOF
LJ-(SR)
Disper.-corr.
Coulomb-(SR)
Coul.-recip.
EOF
    energy_HMX
    energy_ligand

##计算相互作用能
    cp -f  /mmx/gromacs/energy_compute.py ./
    cp -f  /mmx/gromacs/xvg_average.py ./
    ${command1} python energy_compute.py complex.xvg HMX.xvg ligand.xvg -f energy_results.xvg
    result=$($command1 python xvg_average.py energy_results.xvg 5 0 -1 |tail -n 1)
    echo "相互作用能是： $result"
    
# #提取能量
#     echo "12 13\n" |${command1} gmx_mpi energy -f md.edr -o ener_md.xvg 
#     cp -f  /vol8/home/mmx/AI-for-M/data_gmx/scripts/analyze_energy.py $run_path   
#     chmod +x analyze_energy.py
#     result=$($command1 python analyze_energy.py 25)    ##记得根据xvg来进行修改行数，返回能量均值
#     echo $result
#     rm -rf \#*
#     cd $base_path


# ac_energy=$(run_md "ac_${num}" "${ac_path}" |tail -n 1)
# complex_energy=$(run_md "complex_${num}" "${complex_path}" |tail -n 1)
# base_energy=-164672

# delta_energy=$(echo "$complex_energy - $base_energy - $ac_energy" |bc)    #加|bc是进行浮点数运算






