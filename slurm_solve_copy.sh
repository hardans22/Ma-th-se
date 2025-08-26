#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --array=1-2
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=optimum
##SBATCH --nodelist=optimum[04]


# chacune des tâches doit être capable de déterminer qui elle
# est en fonction de son indice dans le array slurm.
# pour cela, on fait toutes les combinaisons possibles
# et si l'indice concorde alors c'est nous et on
# exécute le programme avec les bons arguments.

i=1

for preprocess in false true  
do   
    if [ $SLURM_ARRAY_TASK_ID -eq $i ]
    then
        module load julia
        module load gurobi
        /usr/bin/time -v julia /home/danhar/Documents/Code/__main__copy.jl $preprocess
    fi
    (( i = $i +1 ))
done 

