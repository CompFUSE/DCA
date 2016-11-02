#! /usr/bin/env python3.1

import os
import shutil
import math

print("\n******************************************************** \n")

n        = 16
filename = "trig_ops.h"

file = open(filename, 'w')

file.write("//-*-C++-*-\n\n")
	
file.write("#ifndef COSISNUS_VALUE_H\n")
file.write("#define COSISNUS_VALUE_H\n")
file.write("\n\n")

struct = """template<int i, int n>
struct COSINE_EVAL 
{};\n\n"""

file.write(struct)

struct = """template<int i, int n>
struct SINE_EVAL 
{};\n\n"""

file.write(struct)

struct = """template<>
struct COSINE_EVAL<I_TEMPLATE, N_TEMPLATE>
{
\tstatic double value(){
\t\treturn VALUE;
\t}
};\n\n"""

for i in range(1, n+1):
    for j in range(0, i+1):
        new_struct = struct
        new_struct = new_struct.replace("I_TEMPLATE", str(j))
        new_struct = new_struct.replace("N_TEMPLATE", str(i))
        new_struct = new_struct.replace("VALUE"     , str( math.cos( (2.*math.pi*(j+0.)/(i+0.)) ) ) )
    
        file.write(new_struct)

struct = """template<>
struct SINE_EVAL<I_TEMPLATE, N_TEMPLATE>
{
\tstatic double value(){
\t\treturn VALUE;
\t}
};\n\n"""

for i in range(1, n+1):
    for j in range(0, i+1):
        new_struct = struct
        new_struct = new_struct.replace("I_TEMPLATE", str(j))
        new_struct = new_struct.replace("N_TEMPLATE", str(i))
        new_struct = new_struct.replace("VALUE"     , str(math.sin(2.*math.pi*(j+0.)/(i+0.))))
    
        file.write(new_struct)

file.write("#endif\n\n")

file.close()
