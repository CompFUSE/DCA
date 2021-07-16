#!/bin/bash -l
#SBATCH -J tc_tutorial
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 32
#SBATCH --exclusive
#SBATCH --mem=128g
#SBATCH --gres=gpu:2
#SBATCH -p gpu_p100
#SBATCH -A ccsd
#SBATCH -t 04:00:00


#This will run all the tc calculations
. ../../../build-aux/cades_load_modules.sh
#uncomment the runner apropriate to your system or add one.
RUNNER_CADES="srun -n1 --cpus-per-task=32 --gpus-per-task=1 --accel-bind=gv"

RUN_DCA="${RUNNER_CADES} ../../applications/dca/main_dca"
date

$RUN_DCA ./T=1/input_sp.json
$RUN_DCA ./T=0.75/input_sp.json
$RUN_DCA ./T=0.5/input_sp.json
$RUN_DCA ./T=0.25/input_sp.json
$RUN_DCA ./T=0.125/input_sp.json
$RUN_DCA ./T=0.1/input_sp.json
$RUN_DCA ./T=0.1/input_tp.json
$RUN_DCA ./T=0.09/input_sp.json
$RUN_DCA ./T=0.09/input_tp.json
$RUN_DCA ./T=0.08/input_sp.json
$RUN_DCA ./T=0.08/input_tp.json
$RUN_DCA ./T=0.07/input_sp.json
$RUN_DCA ./T=0.07/input_tp.json

date
