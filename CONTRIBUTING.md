# Contributing

When contributing to CartDG, please adhere to these guidelines. Some are standard (but important) practice
in software development, while others reflect the personal preferences of the maintainer (currently Micaiah)
and are included for the sake of standardization.

## General procedures
Please to **not** commit or push to `master`. The standard workflow should be the following:
1. Create a new branch.
2. Implement changes.
3. Open pull request.

After reviewing the pull request, the maintainer will merge it and delete the branch.
You are encouraged to create as many new branches as you like, for any reason, even if they are
just experiments you don't intend to merge. However, if you are
completely done with a branch, please delete it. You are also encouraged to push frequently (at least
once per work day) -- productivity is maximized when everyone can see the latest version of everything.
Finally, development should be test-driven as much as possible, granted the inherent difficulty of
testing a computational physics solver.

## Pull requests
A pull request will not be merged until it meets the following requirements:
* No merge conflicts.
* Must compile and run in both `Release` and `Debug` mode (the latter with sanitizers).
* No unit test failures.
* No compiler warnings.

## Conventions
Please adhere to the following style conventions whenever possible.
### File formatting
* File names must not contain whitespace or other characters that have meaning in shell commands.
* Source files must not contain the tab character. C++ code should be indented with 2 spaces
  and Python code should be indented with 4 spaces.
* Use blank lines sparingly.
* In functions, loops, etc., the `{` should have a line to itself.

### Variable names
* Should contain words, abbreviations, and/or numbers separated by underscores (not camelcase).
* Names of classes should begin with a capital letter. Otherwise, all letters should be lowercase.
* A name should be a description of the variable it refers to. It
  should not reference mathematical symbols or the sound of words. For example, resist the temptation
  to use `x` to refer to horizontal position. Acceptable alternatives might be `coord0`, or `position`.
* When practical, avoid overly long names. Abbreviations are encouraged, especially in a local scope where there
  are fewer opportunities to confuse similar names.
* Some examples names I hate:
  * `eye` to refer to the identity matrix. Rreferences the sound of the word for the mathematical symbol, two layers removed from any descriptive meaning.
  * `str2num` to mean "string to number". Just call it `str_to_num`. Distinguishing between "2" and "to" is worth 3 extra characters.
  * I have seen names longer than this sentence. Save some for the comments.
## Miscellaneous
* Avoid `using namespace`. The code is more clear when namespaces are referenced explicitly.
* *Everything* should be wrapped in the namespace `cartdg` to avoid creating naming conflicts in NASCART-GT.
  * Macros, which do not respect namespaces, should be prefixed with `CARTDG_`.
* Do not use non-conforming features of GCC.
* The default type for integers should be `int`. For example, don't use `unsigned` for loop indices and array sizes
  unless there is a specific reason to.
