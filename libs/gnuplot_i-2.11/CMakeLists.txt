# Gnuplot C interface

add_library(gnuplot_interface STATIC src/gnuplot_i.c)
target_compile_options(gnuplot_interface PRIVATE "-Wno-implicit-function-declaration")
