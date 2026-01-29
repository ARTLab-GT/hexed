import sys

assert len(sys.argv) == 3
source_dir = sys.argv[1]
target_dir = sys.argv[2]

def get_text(processed):
    text = ""
    for name in sorted(processed.keys()):
        text += f'''\\anchor {name}
<div style="font-size: 150%; padding-bottom: 0.5em;"> __{name}__ </div>
<div style="margin: 0px">
{processed[name]}
</div>
<hr>

'''
    return text

for fname in ["input_parameters", "output_parameters", "macros"]:
    with open(f"{source_dir}/hil/{fname}.hil", "r") as input_file:
        input_text = input_file.read();
    raw = input_text.split('{"')[1:]
    processed = {}
    for r in raw:
        try:
            doc, rest = r.split('"}\n')
        except Exception as e:
            print(r)
            raise e
        if (doc.startswith(r"\variable")):
            name = doc.split(" ")[1].split("\n")[0].strip()
            doc = doc.split("\n")[1:]
            if doc[0].startswith(r"\default"):
                doc[0] = f"__Default:__ `{doc[0][9:]}`\\n"
            else:
                doc = ["__No default value.__\\n"] + doc
            doc = "\n".join(doc)
        else:
            name = rest.split(" =")[0]
            rest = " =".join(rest.split(" =")[1:]).strip()
            if (rest.startswith("{\n")):
                value = f"""<details>
<summary>__Implementation__</summary>
~~~{{.unparsed}}
{rest.split('\n}')[0]}
}}
~~~
</details>
"""
            else:
                value = f"__Default:__ `{rest.split('\n')[0]}`\\n\n"
            doc = value + doc
        processed[name] = doc
    output_text = f"/*! \\page {fname} {fname.replace('_', ' ').title()}\n\n" + get_text(processed) + "*/\n"
    with open(f"{target_dir}/doc/{fname}.dox", "w") as output_file:
        output_file.write(output_text)

with open(f"{source_dir}/libhexed/Case.cpp", "r") as input_file:
    input_text = input_file.read();
    processed = {}
    for command in input_text.split('/*"')[1:]:
        doc, rest = command.split('"*/')
        try:
            name = rest.split('create("')[1].split('"')[0]
        except Exception as e:
            print(rest)
            raise e
        processed[name] = doc.replace("  * ", "")
    output_text = "/*! \\page command_variables Command Variables\n\n" + get_text(processed) + "*/\n"
    with open(f"{target_dir}/doc/command_variables.dox", "w") as output_file:
        output_file.write(output_text)

with open(f"{target_dir}/doc/parameters.dox", "w") as output_file:
    output_file.write(r"""/*! \page parameters Solver Parameters
\subpage input_parameters
\subpage output_parameters
\subpage macros
\subpage command_variables
*/
""")
