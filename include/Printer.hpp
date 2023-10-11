#ifndef HEXED_PRINTER_HPP_
#define HEXED_PRINTER_HPP_

#include <iostream>

namespace hexed
{

//! \brief abstract base class for handling different variations of printing things for user
//! \details basically a slightly higher-level version of output strings
class Printer
{
  public:
  virtual ~Printer() = default;
  virtual void print(std::string) = 0;
};

//! \brief prints to a collection of `std::ostream`s.
class Stream_printer : public Printer
{
  public:
  std::vector<std::ostream*> streams {&std::cout}; //!< \note does not own
  inline void print(std::string message) override {for (auto stream : streams) (*stream) << message << std::flush;}
};

}
#endif
