#!/bin/bash
#
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Running this script remotely triggers a build on Jenkins (project name DCA++-develop).
# The first input argument is the name of the branch or commit to build.
# The script requires curl.
#
# Usage: ./jenkins.sh branch|commit

if [ $# -eq 0 ]; then
    printf "No branch name or commit provided.\nUsage: ./jenkins.sh branch|commit\n"
else
    printf "Triggering build of branch|commit %s.\n" "$1"
    curl https://jenkins.cscs.ch/view/DCA++/job/DCA++-develop-santis/buildWithParameters?token=giro21\&BRANCH=$1    
fi
