#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --array=1-6
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=optimum
##SBATCH --nodelist=optimum[04]

# Génération automatique des instances au lieu de les lister
instances=()
for base in 72 76 83; do
    for nbr_day in 3 ; do
        instances+=("${base}FL_8A_${suffix}D")
    done
done

# chacune des tâches doit être capable de déterminer qui elle
# est en fonction de son indice dans le array slurm.
i=1
for inst_name in "${instances[@]}"; do
    for preprocess in false true; do
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]; then
            module load julia
            module load gurobi
            julia /home/danhar/Documents/Code/__main__.jl $inst_name $preprocess
        fi
        (( i = $i +1 ))
    done
done