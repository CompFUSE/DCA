#! /usr/bin/env python3.1

import commands
import shutil
import os
import sys
import time


print("\n************************************************************************* \n")
print("This will generate dmn_i.h files && include_dmn.h && domain_type_operations")
print("\n************************************************************************* \n")

cmd = "python write_product_domains.py"
os.system(cmd)

cmd = "python write_domain_type_operations.py"
os.system(cmd)

cmd = "python write_include_file.py"
os.system(cmd)




