#ifndef HEXED_CSV_HPP_
#define HEXED_CSV_HPP_

#include <fstream>
#include "Visualizer.hpp"

namespace hexed {

/*! \brief Writes output in Comma-Separated Value format
 * \details File will be open as long as this object exists.
 * If the file already exists, it will be overwritten.
 */
class Csv : public Visualizer {
  int _row;
  int _cols;
  std::ofstream _file;
  public:
  //! \brief Creates a `Csv` with named columns.
  //! \details The resulting CSV file will begin with a line with column names.
  Csv(std::string name, std::vector<std::string> columns);
  //! \brief Creates a `Csv` with unnamed columns.
  //! \details The resulting CSV file will not include a header line, but simply start with the first row of values.
  Csv(std::string name, int n_columns);
  /*! \brief Writes some data to the CSV.
   * \details `data.cols()` must be the same as the number of columns in this file.
   * The number of rows is arbitrary.
   * If you have already written some data to this file with previous calls to `write`
   * (but with this same `Csv` instance), the new data will be appended.
   */
  void write(Array<double> data);
  //! \brief Accepts physical mesh data and writes only the node values. \show_details
  void write_block(Array<double> pos, Array<double> vars) override;
  //! \brief Accepts physical mesh data and writes only the node values. \show_details
  void write_unstruct(Array<Int> elements, Array<double> pos, Array<double> vars) override;
};

}
#endif
