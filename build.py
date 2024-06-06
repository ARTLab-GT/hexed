import subprocess
import sys
subprocess.run(["python3", "python/hexedpy/_build.py"] + sys.argv[1:])
