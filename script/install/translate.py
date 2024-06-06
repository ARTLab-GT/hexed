from math import pi
import re

def translate(in_name, out_name, lang, preamble=""):
    with open(in_name, "r") as in_file:
        text = in_file.read()
    translated = preamble
    if not translated.endswith("\n"):
        translated = translated + "\n"
    if lang == "hil":
        pow_op = "^"
    elif lang == "py":
        pow_op = "**"
    else:
        raise Exception(f"unrecognized language `{lang}`")
    for line in text.split("\n"):
        if re.match(" *const", line): # if this is a const variable declaration
            double = bool(re.search(r"\bdouble\b", line))
            line = re.sub(r" *const +\w+ *", "", line) # remove type
            line = re.sub(";.*", "", line) # remove semicolon an anything after it
            line = re.sub("M_PI", f"{pi:.20f}", line) # replace definition of pi
            # WARNING: the following will only work when both operands are each one word---it won't parse arithmetic expressions
            line = re.sub(r"math::pow\((\w*), ?(\w*)\)", r"\1" + pow_op + r"\2", line) # replace power expressions with operator
            if lang == "hil":
                line = re.sub(r"((?<!\.)\b[0-9]{9,}\b)", r"\1.", line) # `double`ize integer literals that would overflow
                if double: # make sure that things originally intended as doubles stay doubles
                    line += " + 0."
            translated += line + "\n" # add to translated file
    with open(out_name, "w") as out_file:
        out_file.write(translated)

if __name__ == "__main__":
    import sys
    translate(*sys.argv[1:])
