#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=3G
#SBATCH --time=1:10:00
#SBATCH --array=1-100
#SBATCH --output=arrayjob_%A_%a.ou
#SBATCH --partition=cpubase_bycore_b2
#SBATCH --nodelist=c500
#SBATCH --mail-user=harcenagedansou@gmail.com
#SBATCH --mail-type=END


i=1

for Fold in "A_MTN_1" "A_MTN_3" "A_MTN_5" "A_MTN_8" "A_MTN_12"; do
    # Dossier contenant les fichiers JSON
    JSON_DIR="/home/danhar/projects/def-mattgru/danhar/Code/INSTANCES/instances_literature_json/${Fold}"
    for json_file in "$JSON_DIR"/*.json; do
        for sp_method in "MILP" "PD"; do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                module load julia
                module load gurobi
                julia /home/danhar/projects/def-mattgru/danhar/Code/MODEL_FH/BENDERS/__main__benders.jl "$json_file" $sp_method        
                fi
            (( i = $i + 1 ))
        done 
    done
done

