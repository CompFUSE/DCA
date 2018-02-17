# Parallel utils
MESSAGE( "DCA_HAVE_MPI: ${DCA_HAVE_MPI}" )
add_library(parallel_util STATIC get_bounds.cpp ostr_operators.cpp
"$<$<BOOL:${DCA_HAVE_MPI}>:mpi_ostr_operators.cpp>")
