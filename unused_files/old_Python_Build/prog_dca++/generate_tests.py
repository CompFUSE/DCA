import sys
sys.path.append("./python_files/")

import configure

configuration = configure.configure()

#print configuration

configs = [["Todi", "CPU"], 
           ["Todi", "GPU"]]

for l in range(0, len(configs)):

    configuration.write_tests("test_"+configs[l][0]+"_"+configs[l][1], configs[l][0], configs[l][1])

    print configuration
