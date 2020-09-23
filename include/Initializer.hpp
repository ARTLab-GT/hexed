#ifndef INITIALIZER_HPP_
#define INITIALIZER_HPP_

#include <vector>

class Initializer
{
  public:
  Initializer();
  virtual ~Initializer();
  virtual std::vector<double> momentum(std::vector<double> position) = 0;
  virtual std::vector<double> scalar_state(std::vector<double> position) = 0;
};

#endif
