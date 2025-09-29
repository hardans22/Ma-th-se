#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=5:30:00
#SBATCH --array=1-1
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=interac
#SBATCH --nodelist=g4


# Parcourir tous les fichiers JSON

i=1
for json_file in "/home/danhar/projects/def-mattgru/danhar/Code/instances_json/A_MTN_8/249FL_10A_2.json"; do
    for preprocess in false; do
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            module load julia
            module load gurobi
            julia /home/danhar/projects/def-mattgru/danhar/Code/__main__nibi.jl "$ac_file" $preprocess
        fi
        (( i = $i + 1 ))
    done
done
