#!/bin/bash 

#############
#by mmx 20240613  针对complex的top文件产生
#sh fep_gen_top.sh 1  
############


FEP_dir="/mmx/gromacs"
base_dir="${FEP_dir}/NPBA/${1}" 
ready_top_dir="${base_dir}/ready_top"
LIG="LIG_${1}"
HMX="HMX_${1}"
echo $ready_top_dir


mkdir -p "$ready_top_dir"
cd $ready_top_dir
echo "==================="
mv "../$LIG.mol2" $ready_top_dir
mv "../$HMX.mol2" $ready_top_dir
cp -rf  /fs1/home/chengkun_wu/mmx/gromacs/sobtop/* ./


#sobtop程序帮助生成LIG和HMX的.itp、.top、.gro三个文件
gen_sob(){
    local mol2_file="${1}.mol2"   
    if [ ! -f "${1}.top" ];then
        
        ./sobtop   $mol2_file << EOF /dev/null
2
${ready_top_dir}/${1}.gro
1
2
4
${ready_top_dir}/${1}.top
${ready_top_dir}/${1}.itp
0
EOF
    else
        echo "${1}.top文件已存在"
    fi
}

gen_sob "$LIG"
gen_sob "$HMX"


for i in `ls /fs1/home/chengkun_wu/mmx/gromacs/sobtop/`
do
    rm -rf $i 
done

cp -rf $FEP_dir/topol.top ./
#替换16、17行的内容为对应序号分子
sed -i "16s/.*/${HMX}     1/" topol.top
sed -i "17s/.*/${LIG}     1/" topol.top 

#替换每个.gro文件中包含#include行（7，8行）的路径
sed -i  "7s/.*/#include \"${ready_top_dir//\//\\/}\/forcefield.itp\"/" topol.top
sed -i  "8s/.*/#include \"${ready_top_dir//\//\\/}\/$HMX.itp\"/" topol.top
sed -i  "9s/.*/#include \"${ready_top_dir//\//\\/}\/$LIG.itp\"/" topol.top

rm -rf forcefield.itp
touch forcefield.itp
#将HMX和LIG的itp文件的[ atomtypes ]合并至forcefield文件中，并删除1和2.itp文件的[ atomtypes ]部分
sed -n '/\[ atomtypes \]/,/\[ moleculetype \]/{/\[ moleculetype \]/!p;}' $HMX.itp  >forcefield.itp
sed -n '/\[ atomtypes \]/,/\[ moleculetype \]/{/\[ moleculetype \]/!p;}' $LIG.itp  >>forcefield.itp
sed -i '/\[ atomtypes \]/,/\[ moleculetype \]/{/\[ moleculetype \]/!d}' $HMX.itp
sed -i '/\[ atomtypes \]/,/\[ moleculetype \]/{/\[ moleculetype \]/!d}' $LIG.itp
awk '!seen[$0]++' forcefield.itp |sed '/^$/d' >temp.itp
mv temp.itp forcefield.itp



cp -rf $HMX.gro complex.gro
#获取HMX和LIG的.gro文件的原子个数并求和替换complex.gro的原子个数
sum_and_replace(){
    local HMX_gro=$1
    local LIG_gro=$2
    local complex_gro=$3
    #读取HMX和LIG的gro文件的原子个数
    local num_gro1=$(sed -n '2p' $HMX_gro)
    local num_gro2=$(sed -n '2p' $LIG_gro)
    local num_complex=$(($num_gro1 + $num_gro2))
    sed -i "1s/.*/complex/" $complex_gro
    sed -i "2s/.*/${num_complex}/" $complex_gro
    echo "---------复合物原子总个数为${num_complex}-------------"
}

#将complex.gro的最后一行删除，并追加$LIG.gro的原子坐标信息
merge_gro(){
    local complex_gro=$1
    local LIG_gro=$2
    
    sed -i '$d' $complex_gro
    sed -n '3,$p' $LIG_gro >> $complex_gro 
}

sum_and_replace "$HMX.gro" "$LIG.gro" "complex.gro"
merge_gro "complex.gro" "$LIG.gro" 


cd $FEP_dir




