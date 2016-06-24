################################################################################
# Common definitions for building DCA++ programs and tests.
################################################################################

# Libraries
# TODO: Handle git_version and modules libraries nicer.
set(DCA_LIBRARIES "${DCA_EXTERNAL_LIBS};git_version;modules;posix_qmci")
# message("DCA_LIBRARIES: ${DCA_LIBRARIES}")

# Includes
set(DCA_INCLUDES ${DCA_EXTERNAL_INCLUDES})
# message("DCA_INCLUDES: ${DCA_INCLUDES}")
