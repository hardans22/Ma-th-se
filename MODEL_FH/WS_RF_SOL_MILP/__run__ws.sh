#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=35G
#SBATCH --time=05:10:00
#SBATCH --array=1-20
#SBATCH --output=arrayjob_%A_%a.out
#SBATCH --partition=cpubase_bycore_b2
##SBATCH --nodelist=c388
#SBATCH --mail-user=harcenagedansou@gmail.com
#SBATCH --mail-type=END

# Parcourir tous les fichiers JSON

i=1

for Fold in "A_MTN_12" "A_MTN_15"; do
    # Dossier contenant les fichiers JSON
    JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/INSTANCES/instances_literature_json/${Fold}"
    for json_file in "$JSON_DIR"/*.json; do
        for graph_reduc in false; do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                module load julia
                module load gurobi
                julia /home/danhar/projects/def-mattgru/danhar/Code/MODEL_FH/WS_RF_SOL_MILP/__main__ws.jl "$json_file" $graph_reduc
            fi
            ((i = $i + 1))
        done
    done
done 
