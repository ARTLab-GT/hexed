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

all_names = []
page_docs = {"input_parameters": r"""
These are all variables that the solver will read and that will affect its behavior.
Most of them have default values which are initialized before processing your input file,
but a few of them (notably \ref n_dim and \ref reference_length ) are not defined until you assign to them.
You do not need to set all of these---see \ref running for a tutorial on how to write an input file.
""",
"output_parameters": r"""
These are all variables that the solver will assign to.
If you want any information beyond what the solver typically prints to the console,
you can observe these output parameters,
by \ref println "printing" them to the screen or using them in other calculations.
There are also some, such as \ref air_viscosity,
which are simply predefined constants that may be useful in your input files.
""", "solver_macros": r"""
These are builtin \ref Macros which you can execute with the `$` operator.
Some of these, like \ref plot_history and \ref simulate , are intended for you to use directly.
Others, like \ref iterate , are internal macros that are typically called by other builtin macros,
but you can still call them yourself or edit them to change the builtin behavior if you desire.
"""}
for fname in ["input_parameters", "output_parameters", "solver_macros"]:
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
    output_text = f"/*! \\page {fname} {fname.replace('_', ' ').title()}{page_docs[fname]}<hr>\n"
    output_text += get_text(processed) + "*/\n"
    with open(f"{target_dir}/doc/{fname}.dox", "w") as output_file:
        output_file.write(output_text)
    all_names += processed.keys()

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
    output_text = r"""/*! \page command_variables Command Variables
These are \ref heisenberg variables that cause the solver to perform actions when you evaluate them.
Unlike \ref solver_macros, you do not need to use the `$` operator to invoke them;
simply place the name of the variable in its own statement
(which will cause the interpreter to evaluate it but not assign it to anything).
Most of them evaluate to the empty string `{}` so that you can use them like a shell command in an interactive session.
They are typically used internally by the builtin \ref solver_macros, but you can also use them directly.
<hr>
"""
    output_text += get_text(processed) + "*/\n"
    with open(f"{target_dir}/doc/command_variables.dox", "w") as output_file:
        output_file.write(output_text)
    all_names += processed.keys()

with open(f"{target_dir}/doc/parameters.dox", "w") as output_file:
    output_file.write(rf"""/*! \page parameters Solver Parameters
In \ref hil "HIL", as in any other scripting language, you can define as many variables as you like.
However, certain variables names have special meaning.
- \subpage input_parameters are variables that you can change to affect the behavior of the solver.
- \subpage output_parameters are variables that the solver assigns values to in order to communicate information to you.
- \subpage solver_macros are snippets of code that you can execute (or that the solver might execute).
- \subpage command_variables are \ref heisenberg variables
  that cause the solver to perform actions when you evaluate them.

A list of all solver parameters in all four categories is below:
{'\n'.join([r'- \ref ' + name for name in sorted(all_names)])}
*/
""")
