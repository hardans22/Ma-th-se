#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=2:10:00
#SBATCH --array=1-40
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=interac
#SBATCH --nodelist=g4


# Dossier contenant les fichiers JSON
JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/INSTANCES/instances_json/A_MTN_2"

# Parcourir tous les fichiers JSON

i=1

for json_file in "$JSON_DIR"/*.json; do
    if [ $SLURM_ARRAY_TASK_ID -eq $i ]
    then
        module load julia
        module load gurobi
        julia /home/danhar/projects/def-mattgru/danhar/Code/__main__nibi_benders.jl "$json_file"
    fi
    (( i = $i + 1 ))
done

