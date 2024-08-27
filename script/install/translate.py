from math import pi
import re

def translate(in_name, out_name, lang, preamble=""):
    with open(in_name, "r") as in_file:
        text = in_file.read()
    translated = preamble
    if not translated.endswith("\n"):
        translated = translated + "\n"
    namespace = []
    if lang == "hil":
        pow_op = "^"
        def start_class():
            return f"{namespace[-1]} = {{"
        def end_class():
            return "}"
    elif lang == "py":
        pow_op = "**"
        def start_class():
            return f'class {namespace[-1]}:\n    r"""! \\brief stands in for a namespace """'
        def end_class():
            return ""
    else:
        raise Exception(f"unrecognized language `{lang}`")
    for line in text.split("\n"):
        ns_match = re.match(r" *namespace (\w*)", line)
        if ns_match:
            if ns_match.group(1) != "hexed":
                namespace.append(ns_match.group(1))
                translated += "\n" + start_class() + "\n"
        elif re.match(r" *\}", line):
            if namespace:
                translated += end_class() + "\n"
                namespace = namespace[:-1]
        elif re.match(" *const", line): # if this is a const variable declaration
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
            if namespace:
                line = "    "*len(namespace) + line
            translated += line + "\n" # add to translated file
    with open(out_name, "w") as out_file:
        out_file.write(translated)

if __name__ == "__main__":
    import sys
    translate(*sys.argv[1:])
