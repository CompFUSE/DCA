#!/bin/bash

#SBATCH -A cph102
#SBATCH -J DCA_kagome
#SBATCH -o %x-%j.out
#SBATCH -t 00:30:00
#SBATCH -p batch
#SBATCH -q debug
#SBATCH -N 8
#SBATCH --threads-per-core=1
#SBATCH -C nvme

srun --ntasks-per-node 7 -m block:cyclic --gpus-per-task=1 --gpu-bind=closest -n 56 -c 1 ../../applications/dca/main_dca T_1.0/input.sp.in
srun --ntasks-per-node 7 -m block:cyclic --gpus-per-task=1 --gpu-bind=closest -n 56 -c 1 ../../applications/dca/main_dca T_0.4/input.sp.in
srun --ntasks-per-node 7 -m block:cyclic --gpus-per-task=1 --gpu-bind=closest -n 56 -c 1 ../../applications/dca/main_dca T_0.25/input.sp.in
