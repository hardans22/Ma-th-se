#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=7G
#SBATCH --time=5:10:00
#SBATCH --array=1-32
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=interac
#SBATCH --nodelist=g4


# Dossier contenant les fichiers JSON
JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/instances_json/A_MTN_5"

# Parcourir tous les fichiers JSON

i=1

for json_file in "$JSON_DIR"/*.json; do
    for preprocess in false true; do
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            module load julia
            module load gurobi
            julia /home/danhar/projects/def-mattgru/danhar/Code/__main__nibi.jl "$json_file" $preprocess
        fi
        (( i = $i + 1 ))
    done
done

