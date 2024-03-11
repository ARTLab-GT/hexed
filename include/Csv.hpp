#ifndef HEXED_CSV_HPP_
#define HEXED_CSV_HPP_

#include <fstream>
#include "Visualizer.hpp"

namespace hexed
{

class Csv : public Visualizer
{
  int _row;
  int _cols;
  std::ofstream _file;
  public:
  Csv(std::string name, std::vector<std::string> columns);
  void write(Array<double> data);
  void write_block(Array<double> pos, Array<double> vars) override;
  void write_unstruct(Array<int> elements, Array<double> pos, Array<double> vars) override;
};

}
#endif
