#! /usr/bin/env python3.1

import commands
import shutil
import os
import sys
import time

N         = 10
file_name = "domain_type_operations.h"

def print_erase_template(file, n):

    for t1 in range(0,n):
    
        file.write("template <")

        string_i = "class TINDEX "
        for t0 in range(0,n):
            file.write(string_i.replace("INDEX", str(t0)))
            if(t0 == n-1): 
                file.write(">\n")
            else:
                file.write(", ")
        
        file.write("struct Erase<dmn_" + str(n) + "< ")
        
        string_i = "TINDEX"

        for t0 in range(0,n):
            file.write(string_i.replace("INDEX", str(t0)))
            if(t0 == n-1): 
                file.write(">, ")
            else:
                file.write(", ")

        file.write("T"+str(t1)+"> {\n")

        file.write("typedef dmn_" + str(n-1) + "<")

        strs=[]
        for t0 in range(0,n):
            if(t0!=t1):
                strs.append(string_i.replace("INDEX", str(t0)))

        for t0 in range(0,n-1):
            file.write(strs[t0])
            if(t0 == n-2): 
                file.write("> Result;\n")
            else:
                file.write(", ")

        file.write("}; \n\n")

def print_swap_cond_template(file, n):
    
    file.write("template <")

    string_i = "class DINDEX, "
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))

    file.write(" class T1, class T2, int N>\n")
    
    file.write(" struct SWAP_COND<dmn_" + str(n) + "< ")
    
    string_i = "DINDEX"
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n-1): 
            file.write(">, ")
        else:
            file.write(", ")

    file.write(" T1, T2, N> {\n\n")

#     file.write("typedef typename dmn_" + str(n) + "<")

#     string_i = "DINDEX"
#     for t in range(0,n):
#         file.write(string_i.replace("INDEX", str(t)))
#         if(t == n-1): 
#             file.write(">::this_type dmn_type_list;\n\n")
#         else:
#             file.write(", ")

#     file.write("typedef typename TL::Swap<dmn_type_list, T1, T2>::Result swapped_dmn_type_list;\n")

    string_i = "TL::NumberOf<typename DINDEX::this_type, T1>::value"
    string_j = " const static int LENGTH_INDEX = SUM;\n"
    for i_t in range(0,n):

        string_tmp = "N"
        for j_t in range(0,i_t):
            string_tmp = string_tmp + " + " + string_i.replace("INDEX", str(j_t))

        string_tmp = string_j.replace("SUM", string_tmp)
        
        file.write(string_tmp.replace("INDEX", str(i_t)))


    file.write("\ntypedef dmn_" + str(n) + "<\n")

    string_i = "  typename SWAP_COND<DINDEX, T1, T2, LENGTH_INDEX>::Result"
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n-1): 
            file.write("   ")
        else:
            file.write(",\n")

    file.write("> Result;\n")

    file.write("}; \n\n")

def print_swap_first_template(file, n):
    
    file.write("template <")

    string_i = "class DINDEX, "
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))

    file.write(" class T1, class T2>\n")
    
    file.write(" struct SWAP_FIRST<dmn_" + str(n) + "< ")
    
    string_i = "DINDEX"
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n-1): 
            file.write(">, ")
        else:
            file.write(", ")

    file.write(" T1, T2> {\n\n")

    file.write("typedef typename dmn_" + str(n) + "<")

    string_i = "DINDEX"
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n-1): 
            file.write(">::this_type dmn_type_list;\n\n")
        else:
            file.write(", ")

    file.write("typedef typename TL::Swap<dmn_type_list, T1, T2>::Result swapped_dmn_type_list;\n")

    file.write("typedef dmn_" + str(n) + "<\n")

    string_i = " dmn_0<typename TypeAt<swapped_dmn_type_list, INDEX>::Result>"
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n-1): 
            file.write("   ")
        else:
            file.write(",\n")

    file.write("> Result;\n")

    file.write("}; \n\n")

def print_swap_all_template(file, n):
    
    file.write("template <")

    string_i = "class DINDEX, "
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))

    file.write(" class T1, class T2>\n")
    
    file.write(" struct SWAP_ALL<dmn_" + str(n) + "< ")

    string_i = "DINDEX"
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n-1): 
            file.write(">, ")
        else:
            file.write(", ")

    file.write("T1, T2> {       \n")

    file.write("typedef dmn_" + str(n) + "<")

    string_i = "typename SWAP_ALL<DINDEX,T1,T2>::Result"
    for t in range(0,n):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n-1): 
            file.write("   ")
        else:
            file.write(",\n")


    file.write("> Result;\n")

    file.write("}; \n\n")


def print_print_dmn_template(file, n):
    
    file.write("template <")

    string_i = "class DINDEX "
    for t in range(0,n+1):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n): 
            file.write(">\n")
        else:
            file.write(", ")

    file.write("struct printTL<dmn_" + str(n+1) + "<")

    string_i = "DINDEX "
    for t in range(0,n+1):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n): 
            file.write("> > { \n")
        else:
            file.write(", ")


    file.write("static void print() {\n")

    string_i = "printTL<DINDEX>::print(); \n"
    for t in range(0,n+1):
        file.write(string_i.replace("INDEX", str(t)))

    file.write("}\n")

    file.write("};\n\n\n")

def print_print_union_template(file, n):
    
    file.write("template <")

    string_i = "class DINDEX "
    for t in range(0,n+1):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n): 
            file.write(">\n")
        else:
            file.write(", ")

    file.write("struct printTL<union_" + str(n+1) + "<")

    string_i = "DINDEX "
    for t in range(0,n+1):
        file.write(string_i.replace("INDEX", str(t)))
        if(t == n): 
            file.write("> > { \n")
        else:
            file.write(", ")


    file.write("static void print() {\n")

    string_i = "printTL<DINDEX>::print(); \n"
    for t in range(0,n+1):
        file.write(string_i.replace("INDEX", str(t)))

    file.write("}\n")

    file.write("};\n\n\n")





file = open(file_name, 'w')

file.write("//-*-C++-*-\n\n")
file.write("/* \n * \t This is a C++ file generated by write_domain_type_operation.py \n *")

file.write("\n * \t author: Peter Staar  \n */ \n")

file.write("\n#include \"type_list_definitions.h\" \n\n\n")

file.write("#ifndef ")
file.write("DOMAIN_TYPE_OPERATIONS_H_")
file.write("\n")
file.write("#define ")
file.write("DOMAIN_TYPE_OPERATIONS_H_")
file.write("\n\n")

file.write("namespace TL \n{\n")

file.write("\n\n\n")
file.write("/****************************************\n")
file.write("***           ERASE                   ***\n")
file.write("*****************************************/\n")
file.write("\n\n\n")

file.write("template<class T0, class T1> \n")
file.write("struct Erase {}; \n\n")

for l in range(2,N):
    print_erase_template(file, l)

file.write("\n\n\n")
file.write("/****************************************\n")
file.write("***           SWAP-COND               ***\n")
file.write("*****************************************/\n")
file.write("\n\n\n")

file.write("template <class dmn, class T1, class T2, int N> \n")
file.write("struct SWAP_COND {\n")
file.write("  typedef dmn Result;\n")
file.write("};\n\n\n")

file.write("template <class T1,class T2>           \n")
file.write("struct SWAP_COND<dmn_0<T1>, T1, T2, 0> {       \n")
file.write("  typedef dmn_0<T2> Result;            \n")
file.write("};\n\n\n")

file.write("template <class T0, class T1,class T2> \n")
file.write("struct SWAP_COND<dmn_0<T0>, dmn_0<T1>, dmn_0<T2>, 0> {       \n")
file.write("  typedef dmn_0<T0> Result;            \n")
file.write("};\n\n\n")

file.write("template <class T1,class T2>           \n")
file.write("struct SWAP_COND<dmn_0<T1>, dmn_0<T1>, dmn_0<T2>, 0> {       \n")
file.write("  typedef dmn_0<T2> Result;            \n")
file.write("};\n\n\n")

for l in range(1,N):
    print_swap_cond_template(file, l)


file.write("\n\n\n")
file.write("/****************************************\n")
file.write("***           SWAP-FIRST              ***\n")
file.write("*****************************************/\n")
file.write("\n\n\n")

file.write("template <class dmn, class T1,class T2> \n")
file.write("struct SWAP_FIRST\n")
file.write("{\n")
file.write(" typedef typename SWAP_COND<dmn, T1, T2, 0>::Result Result;\n")
file.write("};\n\n\n")

# file.write("template <class dmn, class T1,class T2> \n")
# file.write("struct SWAP_FIRST\n")
# file.write("{};\n\n\n")

# file.write("template <class T0, class T1,class T2> \n")
# file.write("struct SWAP_FIRST<dmn_0<T0>, T1, T2> {       \n")
# file.write("  typedef dmn_0<T0> Result;            \n")
# file.write("};\n\n\n")

# file.write("template <class T1,class T2>           \n")
# file.write("struct SWAP_FIRST<dmn_0<T1>, T1, T2> {       \n")
# file.write("  typedef dmn_0<T2> Result;            \n")
# file.write("};\n\n\n")

# file.write("template <class T0, class T1,class T2> \n")
# file.write("struct SWAP_FIRST<dmn_0<T0>, dmn_0<T1>, dmn_0<T2> > {       \n")
# file.write("  typedef dmn_0<T0> Result;            \n")
# file.write("};\n\n\n")

# file.write("template <class T1,class T2>           \n")
# file.write("struct SWAP_FIRST<dmn_0<T1>, dmn_0<T1>, dmn_0<T2> > {       \n")
# file.write("  typedef dmn_0<T2> Result;            \n")
# file.write("};\n\n\n")

# for l in range(1,N):
#     print_swap_first_template(file, l)


file.write("\n\n\n")
file.write("/****************************************\n")
file.write("***           SWAP-ALL                ***\n")
file.write("*****************************************/\n")
file.write("\n\n\n")

file.write("template <class T0, class T1,class T2> \n")
file.write("struct SWAP_ALL\n")
file.write("{};\n\n\n")

file.write("template <class T0, class T1,class T2> \n")
file.write("struct SWAP_ALL<dmn_0<T0>, T1, T2> {       \n")
file.write("  typedef dmn_0<T0> Result;            \n")
file.write("};\n\n\n")

file.write("template <class T1,class T2>           \n")
file.write("struct SWAP_ALL<dmn_0<T1>, T1, T2> {       \n")
file.write("  typedef dmn_0<T2> Result;            \n")
file.write("};\n\n\n")

file.write("template <class T0, class T1,class T2> \n")
file.write("struct SWAP_ALL<dmn_0<T0>, dmn_0<T1>, dmn_0<T2> > {       \n")
file.write("  typedef dmn_0<T0> Result;            \n")
file.write("};\n\n\n")

file.write("template <class T1,class T2>           \n")
file.write("struct SWAP_ALL<dmn_0<T1>, dmn_0<T1>, dmn_0<T2> > {       \n")
file.write("  typedef dmn_0<T2> Result;            \n")
file.write("};\n\n\n")

for l in range(1,N):
    print_swap_all_template(file, l)
    

file.write("\n\n\n")
file.write("/****************************************\n")
file.write("***           PRINT dmn_i             ***\n")
file.write("*****************************************/\n")
file.write("\n\n\n")

file.write("template <class T0>                       \n")
file.write("struct printTL<dmn_0<T0> > {              \n")
file.write("\tstatic void print() {                   \n")
file.write("\t\tstd::cout <<  \"\\t\" << __PRETTY_FUNCTION__ << \"\\n\"; \n");
file.write("\t}                                       \n")
file.write("};                                    \n\n\n")

for l in range(0,N):
    print_print_dmn_template(file, l)


file.write("\n\n\n")
file.write("/****************************************\n")
file.write("***           PRINT union_i           ***\n")
file.write("*****************************************/\n")
file.write("\n\n\n")

for l in range(0,N):
    print_print_union_template(file, l)

file.write("\n}\n")
file.write("\n\n#endif\n") 			
	
file.close()
