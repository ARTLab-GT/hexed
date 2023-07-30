#include <Case.hpp>

namespace hexed
{

Case::Case(std::string input_file)
{
  inter.exec(format_str(1000, "$read \"%s/include/Case.hil\"", config::root_dir));
  inter.exec(format_str(1000, "$read \"%s\"", input_file.c_str()));
}

}
