#include <hexed/config.hpp>
#include <hexed/Case.hpp>

int main(int argc, char* argv[])
{
  #if HEXED_THREADED
  omp_set_num_threads(omp_get_max_threads() - 2);
  #endif
  if (argc == 1) hexed::Case c;
  else if (argc == 2) hexed::Case c(argv[1]);
  else throw std::runtime_error("`hexecute` accepts exactly 0 or 1 arguments");
}
