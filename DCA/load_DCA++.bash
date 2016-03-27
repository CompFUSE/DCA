#!/bin/bash
#
# Sources the easybuild setup script specifiying /project/s299/easybuild/daint
# as the project directory and afterwards loads the DCA++ module created with
# easybuild.

#if the arguments are provided loas eb into from $1 and build into $2
#else revert to standard for daint.
if [ "$#" -eq 3 ]; then
EB_SCRIPT=$1
EB_FOLDER=$2
else 
echo "Using default folders for daint"
EB_SCRIPT=/apps/common/easybuild/setup.sh
EB_FOLDER=/project/s299/easybuild/daint 
fi
echo "Loading easybild from $EB_SCRIPT into $EB_FOLDER"

source $EB_SCRIPT  $EB_FOLDER
DCA_SOURCE="$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "running scripts in ${DCA_SOURCE}/easybuild"
export EASYBUILD_ROBOT_PATHS=$DCA_SOURCE/easybuild:$EASYBUILD_ROBOT_PATHS
eb $DCA_SOURCE/easybuild/DCA++-CMake-CrayGNU-2015.11.eb -r
echo ""
echo "Loading DCA++ ..."
module load DCA++/1.0-CrayGNU-2015.11
