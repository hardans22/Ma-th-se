#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --mem=250G
#SBATCH --time=41:00:00
#SBATCH --array=1-8
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=optimum
##SBATCH --nodelist=optimum[04]


# chacune des tâches doit être capable de déterminer qui elle
# est en fonction de son indice dans le array slurm.
# pour cela, on fait toutes les combinaisons possibles
# et si l'indice concorde alors c'est nous et on
# exécute le programme avec les bons arguments.

i=1
for inst_name in 89FL_7A_3D 171FL_16A_3D 186FL_7A_7D 370FL_15A_7D
do
    for preprocess in false true  
    do   
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            module load julia
            module load gurobi
            /usr/bin/time -v julia /home/danhar/Documents/Code/__main__.jl $inst_name $preprocess
        fi
        (( i = $i +1 ))
    done 
done
 
