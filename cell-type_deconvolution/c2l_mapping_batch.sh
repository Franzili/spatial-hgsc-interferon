#!/bin/bash
#SBATCH --job-name=c2l_mapping
#SBATCH --account=project_xxxxxxx
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --gres=gpu:v100:4,nvme:8
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load pytorch
# Activate the project Python environment for your cluster, e.g.:
#   source "${ENV_PATH}/activate"
source "${ENV_PATH:?set ENV_PATH to the bin directory of your Python environment}/activate"

srun python c2l_spatial_mapping_zheng.py
