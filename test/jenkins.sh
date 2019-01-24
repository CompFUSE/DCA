#!/bin/bash
#
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# This script remotely triggers a build on Jenkins (project name DCA++-develop).
# Running the script requires a Jenkins account and the script will ask for username and password.
# Dependencies: curl.
#
# Usage: ./jenkins.sh

printf "Branch name or commit hash: "
read BRANCH

printf "Username: "
read USERNAME

stty -echo
printf "Password: "
read PASSWORD
stty echo
printf "\n"

printf "Triggering build of branch|commit $BRANCH.\n"
curl -u $USERNAME:$PASSWORD https://ci.cscs.ch:7000/job/s299/job/DCA++-develop/buildWithParameters?token=xjs2ks9m\&BRANCH=$BRANCH
