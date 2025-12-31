#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=02:10:00
#SBATCH --array=1-160
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=interac
#SBATCH --nodelist=c[170]
#SBATCH --mail-user=harcenagedansou@gmail.com
#SBATCH --mail-type=END

# Dossier contenant les fichiers JSON
JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/INSTANCES/instances_json/A_MTN_2"

# Parcourir tous les fichiers JSON

i=1

for json_file in "$JSON_DIR"/*.json; do
    for preprocess in false; do
        for option in MIN MAX MEDIAN MEAN; do 
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                module load julia
                module load gurobi
                julia /home/danhar/projects/def-mattgru/danhar/Code/SENSITIVE_EXP/main_exp_nibi.jl "$json_file" $preprocess $tk_option $option 
            fi
            ((i = $i + 1))
        done 
    done
done
