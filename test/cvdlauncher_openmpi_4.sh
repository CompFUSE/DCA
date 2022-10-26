#!/bin/bash
# simple wrapper setting CUDA_VISIBLE_DEVICES for jobs

export CUDA_VISIBLE_DEVICES=$((OMPI_COMM_WORLD_LOCAL_RANK/4))
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"

$cmd $*
