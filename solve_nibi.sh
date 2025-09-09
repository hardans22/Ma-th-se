#!/bin/bash
#SBATCH --cpus-per-task=10
##SBATCH --mem=10G
#SBATCH --time=6:00:00
#SBATCH --array=1-128
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition= interac
##SBATCH --nodelist=optimum[04]


# Dossier contenant les fichiers JSON
JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/instances_json"

# Parcourir tous les fichiers JSON

i=1
for week_folder in "$JSON_DIR"/M01_W*; do
    for ms_folder in "$week_folder"/MS_*; do
        for json_file in "$ms_folder"/*.json; do
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
    done
done

