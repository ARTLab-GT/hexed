#include <Case.hpp>
#include <iostream>

namespace hexed
{

Case::Case(std::string input_file)
{
  inter.exec(format_str(1000, "$read \"%s\"", input_file.c_str()));
}

}
