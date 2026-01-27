import sys

assert len(sys.argv) == 3
source = sys.argv[1]
target = sys.argv[2]

with open(source, "r") as input_file:
    input_text = input_file.read();
raw = input_text.split('{"')[1:]
processed = {}
for r in raw:
    doc, rest = r.split('"}\n')
    name, rest = rest.split(" =")[:2]
    doc = f"__Default:__ `{rest.split('\n')[0].strip()}`\\n\n{doc}"
    processed[name] = doc
output_text = "/*! \\page parameters Solver Parameters\n"
for name in sorted(processed.keys()):
    output_text += f'''\n<div style="font-size: 150%; padding-bottom: 0.5em;"> __{name}__ </div>
<div style="margin: 0px">
\\anchor {name}
{processed[name]}
</div>
<hr>\n'''
output_text += "*/\n"
with open(target, "w") as output_file:
    output_file.write(output_text)
