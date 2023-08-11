#include <hexed/Case.hpp>

int main(int argc, char* argv[])
{
  if (argc != 2) throw std::runtime_error("`hexed_exec` takes exactly one argument, namely the path to the input file");
  hexed::Case c(argv[1]);
}
