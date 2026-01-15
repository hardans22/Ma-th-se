#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --array=1-10
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=interac
#SBATCH --nodelist=c[170]
#SBATCH --mail-user=harcenagedansou@gmail.com
#SBATCH --mail-type=END

# Dossier contenant les fichiers JSON
JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/INSTANCES/instances_literature_json/A_MTN_1"

# Parcourir tous les fichiers JSON

i=1

for json_file in "$JSON_DIR"/*.json; do
    for graph_reduc in false; do
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            module load julia
            module load gurobi
            julia /home/danhar/projects/def-mattgru/danhar/Code/NORMAL/main_RF_nibi.jl "$json_file" $graph_reduc
        fi
        ((i = $i + 1))
    done
done
