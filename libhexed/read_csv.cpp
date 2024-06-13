#include <fstream>
#include <read_csv.hpp>

namespace hexed
{

Eigen::Matrix<double, dyn, dyn, Eigen::RowMajor> read_csv(std::string file_name)
{
  std::vector<double> data; // use a vector for accumulating the data as i suspect it may handle reallocation more efficiently
  std::ifstream file(file_name);
  HEXED_ASSERT(file.good(), format_str(1000, "failed to open file `%s`", file_name.c_str()));
  std::string line;
  std::getline(file, line);
  int cols = std::count(line.begin(), line.end(), ',') + 1;
  HEXED_ASSERT(cols > 0, "nonpositive number of columns");
  int rows = 0;
  while (!line.empty()) {
    ++rows;
    auto begin = line.begin();
    for (int i_col = 0; i_col < cols; ++i_col) {
      auto end = std::find(begin, line.end(), ',');
      std::string entry(begin, end);
      try {
        data.push_back(std::stod(entry));
      } catch (...) {
        HEXED_ASSERT(false, format_str(1000, "could not convert `%s` to `double`", entry.c_str()))
      }
      begin = end + 1;
    }
    std::getline(file, line);
  }
  Eigen::Matrix<double, dyn, dyn, Eigen::RowMajor> mat(rows, cols);
  for (unsigned i = 0; i < data.size(); ++i) mat(i) = data[i];
  return mat;
}

}
