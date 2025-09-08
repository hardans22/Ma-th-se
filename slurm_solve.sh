#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --mem=200G
#SBATCH --time=85:00:00
#SBATCH --array=1-8
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=optimumlong
##SBATCH --nodelist=optimum[04]

# Génération automatique des instances au lieu de les lister

# chacune des tâches doit être capable de déterminer qui elle
# est en fonction de son indice dans le array slurm.
i=1
for nbr_ac in 5 10; do
    for nbr_day in 7 14; do
        for preprocess in false true; do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]; then
                module load julia
                module load gurobi
                /usr/bin/time -v julia /home/danhar/Documents/Code/__main__slurm.jl $nbr_ac $nbr_day $preprocess
            fi
            (( i = $i +1 ))
        done
    done
done