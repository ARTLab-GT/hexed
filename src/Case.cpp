#include <filesystem>
#include <Case.hpp>

namespace hexed
{

Solver& Case::_solver()
{
  HEXED_ASSERT(_solver_ptr, "`Solver` object does not exist");
  return *_solver_ptr;
}

Case::Case(std::string input_file)
{
  // create custom Heisenberg variables
  _inter.variables->create<int>("create_solver", new Namespace::Heisenberg<int>([this]() {
    auto n_dim = _inter.variables->lookup<int>("n_dim");
    HEXED_ASSERT(n_dim && n_dim.value() > 0 && n_dim.value() <= 3,
                 "`n_dim` must be defined as an integer in [1, 3]");
    auto row_size = _inter.variables->lookup<int>("row_size");
    HEXED_ASSERT(row_size.value() >= 2 && row_size.value() <= config::max_row_size,
                 format_str(300, "`row_size` must be between 2 and %i", config::max_row_size));
    _solver_ptr.reset(new Solver(*n_dim, *row_size, 1.));
    return 0;
  }));
  _inter.variables->create<int>("visualize", new Namespace::Heisenberg<int>([this]() {
    std::filesystem::path file_name(_inter.variables->lookup<std::string>("vis_file_name").value());
    if (!std::filesystem::exists(file_name.parent_path())) {
      std::filesystem::create_directory(file_name.parent_path());
    }
    _solver().visualize_field_tecplot(file_name.string());
    return 0;
  }));
  // load HIL code for the Case _interface
  _inter.exec(format_str(1000, "$read \"%s/include/Case.hil\"", config::root_dir));
  // execute input file
  _inter.exec(format_str(1000, "$read \"%s\"", input_file.c_str()));
}

}
