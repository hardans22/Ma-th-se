#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=02:10:00
#SBATCH --array=1-40
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=interac
#SBATCH --nodelist=g4
#SBATCH --mail-user=harcenagedansou@gmail.com
##SBATCH --mail-type=END

# Dossier contenant les fichiers JSON
#JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/INSTANCES/instances_json/A_MTN_2"

# Parcourir tous les fichiers JSON

i=1
for Fold in "A_MTN_1" "A_MTN_3" "A_MTN_5" "A_MTN_8"; do
    # Dossier contenant les fichiers JSON
    JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/INSTANCES/instances_literature_json/${Fold}"
    for json_file in "$JSON_DIR"/*.json; do
        for preprocess in false; do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                module load julia
                module load gurobi
                julia /home/danhar/projects/def-mattgru/danhar/Code/MODEL_FH/__main__milp_nibi.jl "$json_file" $preprocess
            fi
            ((i = $i + 1))
        done
    done
done