################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Common definitions for building DCA++ applications and tests.

# Libraries
set(DCA_LIBRARIES ${DCA_EXTERNAL_LIBS} git_version modules posix_qmci gaussian_quadrature random)
# message("DCA_LIBRARIES = ${DCA_LIBRARIES}")

# Includes
set(DCA_INCLUDES ${DCA_EXTERNAL_INCLUDES})
# message("DCA_INCLUDES = ${DCA_INCLUDES}")
