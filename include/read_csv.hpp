#ifndef HEXED_READ_CSV_HPP_
#define HEXED_READ_CSV_HPP_

#include "math.hpp"

namespace hexed
{

/*! \brief Reads a CSV file to a matrix.
 * \details Reads a table of data from an ASCII file in Comma Separated Value format.
 * Every row must have the same number of columns.
 * Each entry must contain exactly one numeric literal,
 * preceded and/or followed by zero or more spaces and/or tabs, and nothing else.
 * Thus trailing commas are not allowed, as this would be construed as an empty column.
 * File extension does not matter.
 */
Eigen::Matrix<double, dyn, dyn, Eigen::RowMajor> read_csv(std::string file_name);

}
#endif
