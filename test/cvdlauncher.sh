#!/bin/bash
# simple wrapper setting CUDA_VISIBLE_DEVICES for jobs

let ngpus=6;
export cmd=$1
if [ ! -z $CUDA_VISIBLE_DEVICES ]
then
   let gpu=$CUDA_VISIBLE_DEVICES
   let mydevice=$gpu
   let gpu_start=$gpu+1
  for ((g=$gpu_start; g<$ngpus; g++))
   do
     mydevice=${mydevice},$g
   done

   for ((g=0; g<$gpu; g++))
   do
    mydevice=${mydevice},$g
   done
   export CUDA_VISIBLE_DEVICES=$mydevice
fi
shift
$cmd $*
