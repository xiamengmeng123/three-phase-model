#!/bin/bash

# Set some environment variables 
partition=$1
base_dir=$3
cd $base_dir
FREE_ENERGY=`pwd`
echo "Free energy home directory set to $FREE_ENERGY"
MDP=$FREE_ENERGY/MDP
echo ".mdp files are stored in $MDP"



INIGRO=complex-equil.gro
TOPOL=complex-res.top

#MDRUN="gmx_mpi mdrun -ntmpi 1 -ntomp 6" #11 ns/day

# Change to the location of your GROMACS-2018 installation
#NPROC=1

#Pfifo="/tmp/$$.fifo"
#mkfifo $Pfifo
#exec 6<>$Pfifo
#rm -f $Pfifo

#for ((i=1;i<=$NPROC;i++)); do
#        echo
#done >&6


DIR_ROOT=`pwd`
MOLS=`ls`

export OMP_NUM_THREADS=1

#    read -u6
#    GPU_ID=`getgpu.sh|head -n 1`
#    {
#    gpu_id_local=$GPU_ID
    #MDRUN="gmx_mpi mdrun -ntmpi 3 -ntomp 2 -gpu_id $gpu_id_local" 
    #MDMIN="gmx_mpi mdrun -ntmpi 1 -ntomp 7 -nb gpu -gpu_id $gpu_id_local"
    #MDRUN="gmx_mpi mdrun -ntmpi 1 -ntomp 7 -nb gpu -bonded gpu -gpu_id $gpu_id_local"
MDMIN="yhrun -N 1  -p $partition"
MDRUN="yhrun -N 2 -n 72 -p $partition  gmx_mpi mdrun"
LAMBDA=$2

# A new directory will be created for each value of lambda and
# at each step in the workflow for maximum organization.

mkdir Lambda_$LAMBDA
cd Lambda_$LAMBDA
traj_file="Production_MD/md$LAMBDA.xvg"
if [ -f "$traj_file" ] && tail -n 1  "$traj_file" |grep  "^4000\.0000"; then
    echo "md$LAMBDA.xvg轨迹文件已存在"
else
    echo "md$LAMBDA.xvg文件不完整，重新run"
    rm -rf step*.pdb
    ##############################
    # ENERGY MINIMIZATION STEEP  #
    ##############################
    echo "Starting minimization for lambda = $LAMBDA..." 

    mkdir EM
    cd EM

    # Iterative calls to grompp and mdrun to run the simulations

    $MDMIN gmx_mpi grompp -f $MDP/em_steep_$LAMBDA.mdp -c $FREE_ENERGY/$INIGRO -p $FREE_ENERGY/$TOPOL -o min$LAMBDA.tpr -maxwarn 100

    $MDMIN gmx_mpi mdrun  -v -deffnm min$LAMBDA

    rm -rf \#*

    #####################
    # NVT EQUILIBRATION #
    #####################
    echo "Starting constant volume equilibration..."

    cd ../
    mkdir NVT
    cd NVT

    $MDMIN gmx_mpi grompp -f $MDP/nvt_$LAMBDA.mdp -c ../EM/min$LAMBDA.gro -r ../EM/min$LAMBDA.gro -p $FREE_ENERGY/$TOPOL  -o nvt$LAMBDA.tpr -maxwarn 100

    $MDRUN -v -deffnm nvt$LAMBDA
    rm -rf \#*
    echo "Constant volume equilibration complete."

    #####################
    # # NPT EQUILIBRATION #
    # #####################
    # echo "Starting constant pressure equilibration..."

    # cd ../
    # mkdir -p NPT
    # cd NPT

    # gmx grompp -f $MDP/npt_$LAMBDA.mdp -c ../NVT/nvt$LAMBDA.gro -r ../NVT/nvt$LAMBDA.gro -p $FREE_ENERGY/$TOPOL -n $FREE_ENERGY/index.ndx -t ../NVT/nvt$LAMBDA.cpt -o npt$LAMBDA.tpr -maxwarn 100

    # $MDRUN -deffnm npt$LAMBDA

    # echo "Constant pressure equilibration complete."

    #################
    # PRODUCTION MD #
    #################
    echo "Starting production MD simulation..."

    cd ../
    mkdir Production_MD
    cd Production_MD

    $MDMIN gmx_mpi grompp -f $MDP/md_$LAMBDA.mdp -c ../NVT/nvt$LAMBDA.gro -p $FREE_ENERGY/$TOPOL  -t ../NVT/nvt$LAMBDA.cpt -o md$LAMBDA.tpr -maxwarn 100

    $MDRUN -v -deffnm md$LAMBDA
    rm -rf \#*
    echo "Production MD complete."

    # End
    echo "Ending. Job completed for lambda = $LAMBDA"

    cd $FREE_ENERGY
fi
#    sleep 1
#    echo >&6
#    }&
#    sleep 10


#wait

#exec 6>&-

