import sys

assert len(sys.argv) == 3
source = sys.argv[1]
target = sys.argv[2]

with open(source, "r") as input_file:
    input_text = input_file.read();
output_text = r"/*! \page parameters Solver Parameters"
output_text += "*/\n"
with open(target, "w") as output_file:
    output_file.write(output_text)
