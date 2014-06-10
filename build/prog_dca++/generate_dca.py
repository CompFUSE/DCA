import sys
sys.path.append("./python_files/")

import configure

configuration = configure.configure()

print configuration

configuration.read()

print configuration

configuration.write_all()
