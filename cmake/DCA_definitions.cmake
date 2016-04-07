################################################################################
# Common definitions for building DCA++ programs and tests.
################################################################################

# Libraries
# FIXME: Handle gitVersion and modules library nicer.
set(DCA_LIBRARIES "${DCA_EXTERNAL_LIBS};gitVersion;modules")
# message("DCA_LIBRARIES: ${DCA_LIBRARIES}")

# Includes
set(DCA_INCLUDES
  ${PROJECT_SOURCE_DIR}/src
   )
list(APPEND DCA_INCLUDES "${PROJECT_SOURCE_DIR}/gitVersion")  # Directory of gitVersion.hpp.
list(APPEND DCA_INCLUDES "${PROJECT_SOURCE_DIR}/modules")  # Directory of modules.hpp.
list(APPEND DCA_INCLUDES "${DCA_EXTERNAL_INCLUDES}")
# message("DCA_INCLUDES: ${DCA_INCLUDES}")
