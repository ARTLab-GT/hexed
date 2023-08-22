#include <hexed/Case.hpp>

int main(int argc, char* argv[])
{
  if (argc == 1) hexed::Case c;
  else if (argc == 2) hexed::Case c(argv[1]);
  else throw std::runtime_error("`hexecute` accepts exactly 0 or 1 arguments");
}
