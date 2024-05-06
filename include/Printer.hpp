#ifndef HEXED_PRINTER_HPP_
#define HEXED_PRINTER_HPP_

#include <iostream>
#include "utils.hpp"

namespace hexed
{

//! \brief abstract base class for handling different variations of printing things for user
//! \details basically a slightly higher-level version of output strings
class Printer
{
  public:
  virtual ~Printer() = default;
  virtual void operator()(std::string, bool emph = false) = 0;
};

class Compound_printer : public Printer
{
  public:
  std::vector<std::shared_ptr<Printer>> printers;
  inline void operator()(std::string message, bool emph = false) override {for (auto& printer : printers) (*printer)(message, emph);}
};

//! \brief prints to a `std::ostream`.
class Stream_printer : public Printer
{
  std::ostream& _stream;
  std::string _format_code;
  std::string _reset_code;
  public:
  enum format_type {
    unspecified_type = -1,
    bold = 1,
    dim = 2,
    underline = 4,
    blink = 5,
    reverse = 7,
    hidden = 8,
  };
  enum format_color {
    unspecified_color = -1,
    default_color = 9,
    black = 0,
    red = 1,
    green = 2,
    yellow = 3,
    blue = 4,
    magenta = 5,
    cyan = 6,
    gray = 7,
  };
  struct Format {
    format_type type;
    format_color color;
    bool light;
    bool background;
  };
  Stream_printer(std::ostream& stream = std::cout, Format emph_format = {
    .type = unspecified_type,
    .color = unspecified_color,
    .light = false,
    .background = false})
  : _stream{stream}, _format_code(""), _reset_code("")
  {
    if (!(emph_format.type == unspecified_type && emph_format.color == unspecified_color)) {
      _format_code = format_str(100, "\x1b[%i;%i%im", emph_format.type, 3 + 6*emph_format.light + emph_format.background, emph_format.color);
      _reset_code = "\x1b[0m";
    }
  }
  inline void operator()(std::string message, bool emph = false) override
  {
    _stream << (emph ? _format_code : "") << message << (emph ? _reset_code : "") << std::flush;
  }
};

struct Printer_set
{
  Compound_printer info;
  Compound_printer warn;
  Compound_printer error;
  Printer_set()
  {
    info.printers.emplace_back(std::make_shared<Stream_printer>());
    warn.printers.emplace_back(std::make_shared<Stream_printer>(std::cerr, Stream_printer::Format{.color = Stream_printer::yellow, .light = true}));
    error.printers.emplace_back(std::make_shared<Stream_printer>(std::cerr, Stream_printer::Format{.type = Stream_printer::bold,
                                                                                                   .color = Stream_printer::red}));
  }
};


}
#endif
