#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=01:10:00
#SBATCH --array=1-20
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=cpubase_bycore_b2
#SBATCH --nodelist=c388
#SBATCH --mail-user=harcenagedansou@gmail.com
#SBATCH --mail-type=BEGIN,END

# Parcourir tous les fichiers JSON

i=1

for Fold in "A_MTN_8" "A_MTN_12"; do
    # Dossier contenant les fichiers JSON
    JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/INSTANCES/instances_new_json/${Fold}"
    for json_file in "$JSON_DIR"/*.json; do
        for graph_reduc in false; do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                module load julia
                module load gurobi
                julia /home/danhar/projects/def-mattgru/danhar/Code/MODEL_FH/RF_FO/__main__rffo.jl "$json_file" $graph_reduc
            fi
            ((i = $i + 1))
        done
    done
done 
