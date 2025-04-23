#include <Csv.hpp>
#include <utils.hpp>

namespace hexed {

Csv::Csv(std::string name, int n_columns) : _row{0}, _cols{n_columns}, _file(name + ".csv") {}

Csv::Csv(std::string name, std::vector<std::string> columns) : Csv(name, columns.size()) {
  std::string text;
  for (auto col : columns) text += col + ",";
  text.pop_back();
  _file << text << "\n";
}

void Csv::write(Array<double> data) {
  HEXED_ASSERT(data.order() == 2, "CSV data must be 2D");
  auto shape = data.shape();
  HEXED_ASSERT(shape[1] == _cols, "CSV data has wrong number of columns");
  for (int i_row = 0; i_row < shape[0]; ++i_row) {
    std::string text;
    for (int col = 0; col < shape[1]; ++col) text += format_str(30, "%.20e,", data(i_row)[col]);
    text.pop_back();
    _file << text << "\n";
    ++_row;
  }
}

void Csv::write_block(Array<double> pos, Array<double> vars) {
  #if 0
  HEXED_ASSERT(pos.order() > 1, "input arrays have wrong order");
  pos.reshape({same, whatever});
  vars.reshape({same, whatever});
  HEXED_ASSERT(pos.shape()[1] = vars.shape()[1], "`pos` and `vars` must have compatible shape")
  HEXED_ASSERT(pos.shape()[0] + vars.shape()[0] == _cols,
               "total number of position and state variables must equal number of columns")
  int rows = pos(0).size();
  for (int i_row = 0; i_row < rows; ++i_row) {
    std::string text;
    for (int col = 0; col <  pos.shape()[0]; ++col) text += format_str(30, "%.20e,", pos(col)[i_row]);
    for (int col = 0; col < vars.shape()[0]; ++col) text += format_str(30, "%.20e,", vars(col)[i_row]);
    text.pop_back();
    _file << text << "\n";
    ++_row;
  }
  #endif
}

void Csv::write_unstruct(Array<Int> elements, Array<double> pos, Array<double> vars) {
  write_block(pos, vars);
}

}
