#include <hexed/Interpreter.hpp>

int main(int argc, char* argv [])
{
  if (argc < 1 || argc > 2) throw std::runtime_error("`hil` takes 0 or 1 arguments");
  hexed::Interpreter inter;
  if (argc == 2) inter.exec("$read " + std::string(argv[1]));
  else inter.exec("$repl");
  return 0;
}
