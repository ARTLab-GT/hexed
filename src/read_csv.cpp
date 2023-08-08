#include <fstream>
#include <read_csv.hpp>

namespace hexed
{

Eigen::Matrix<double, dyn, dyn, Eigen::RowMajor> read_csv(std::string file_name)
{
  return Mat<>::Zero(1, 1);
}

}
