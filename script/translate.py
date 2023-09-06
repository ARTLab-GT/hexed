import sys
import re
from math import pi

assert len(sys.argv) == 2, "translate.py requires exactly 1 argument: name of file to process"
file_name = sys.argv[1]
with open(file_name, "r") as in_file:
    text = in_file.read()

hil_text = ""
for line in text.split("\n"):
    if re.match(" *const", line): # if this is a const variable declaration
        double = bool(re.search(r"\bdouble\b", line))
        line = re.sub(" *const +\w+ *", "", line) # remove type
        line = re.sub(";.*", "", line) # remove semicolon an anything after it
        line = re.sub("M_PI", f"{pi:.20f}", line) # replace definition of pi
        # WARNING: the following will only work when both operands are each one word -- it won't parse arithmetic expressions
        line = re.sub("math::pow\((\w*), ?(\w*)\)", r"\1^\2", line) # replace power expressions with `^` operator
        line = re.sub(r"((?<!\.)\b[0-9]{9,}\b)", r"\1.", line) # `double`ize integer literals that would overflow
        if double: # make sure that things originally intended as doubles stay doubles
            line += " + 0."
        hil_text += line + "\n" # add to HIL script

with open(file_name.split("/")[-1].replace(".hpp", ".hil"), "w") as out_file:
    out_file.write(hil_text)
