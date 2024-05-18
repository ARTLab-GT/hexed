#ifndef HEXED_PRINTER_HPP_
#define HEXED_PRINTER_HPP_

#include <iostream>
#include <vector>
#include <memory>

namespace hexed
{

//! \brief abstract base class for handling different variations of printing things for user
//! \details basically a slightly higher-level version of output streams
class Printer
{
  public:
  virtual ~Printer() = default;
  //! \param message what you want to print
  //! \param emph if `true`, some form of emphasis formatting will be used, if possible
  virtual void operator()(std::string message, bool emph = false) = 0;
};

//! \brief `Printer` which is just a combination of other printers.
//! \details Anything you tell it to print will be forwarded to all of `Compound_printer::printers`
class Compound_printer : public Printer
{
  public:
  std::vector<std::shared_ptr<Printer>> printers; //!< \brief list of printers---feel free to modify
  inline void operator()(std::string message, bool emph = false) override {for (auto& printer : printers) (*printer)(message, emph);}
};

//! \brief prints to a `std::ostream`.
class Stream_printer : public Printer
{
  std::ostream& _stream;
  std::string _format_code;
  std::string _reset_code;

  public:
  //! \brief specifies the typeface (i.e. font) to print in
  enum format_type {
    unspecified_type = -1, //!< \brief indicates that the type is not specified
    default_type = 0,
    bold = 1,
    dim = 2,
    underline = 4,
    blink = 5,
    reverse = 7, //!< \brief indicates that the foreground and background colors are swapped
    hidden = 8,
  };

  /*! \brief specifies the color to print in
   * \note `unspecified_color` is not the same as `default_color`.
   * The latter will explicitly send an [ASCII escape code](https://en.wikipedia.org/wiki/ANSI_escape_code#Colors) which sets the color to the default.
   * If the output stream was already printing in some non-default color, this will reset it to the default.
   * It also requires the output stream to support colors.
   * On the other hand, `unspecified_color` indicates that no escape code will be sent to specify the color
   */
  enum format_color {
    unspecified_color = -1, //!< \brief indicates that the text color is not specified
    default_color = 9, //!< \brief explicitly indicates that the default color should be used
    black = 0,
    red = 1,
    green = 2,
    yellow = 3,
    blue = 4,
    magenta = 5,
    cyan = 6,
    gray = 7,
  };

  //! \brief for passing text formatting parameters to `Stream_printer::Stream_printer`
  struct Format {
    format_type type = unspecified_type;
    format_color color = unspecified_color;
    bool light = false; //!< \brief if `true`, make the colors lighter than they otherwise would be
    bool background = false; //!< \brief if `true`, any colors specified refer to the text background rather than the text itself
  };
  static const Format default_format;

  /*! \brief constructs a `Stream_printer`
   * \param stream stream to print to
   * \param emph_format Formatting to apply to emphasized text.
   *   This requires `stream` to support [ASCII formatting codes](https://en.wikipedia.org/wiki/ANSI_escape_code#Colors).
   *   `cout` and `cerr` support these codes, but ASCII files do not.
   */
  Stream_printer(std::ostream& stream = std::cout, Format emph_format = default_format);

  void operator()(std::string, bool emph = false) override;
};

//! \brief A complete set of printers that support messages with various purposes
struct Printer_set
{
  Compound_printer info; //!< \brief general information
  Compound_printer warn; //!< \brief (non-fatal) warnings
  Compound_printer error; //!< \brief (fatal) errors
  /*! \brief Constructs a `Printer_set` with printers initialized to standard output/error
   * \details You can later modify the individual printers.
   * Specifically, the initial values are:
   * - `info` is `std::cout` without any special formatting
   * - `warn` is `std::cerr` with light yellow emphasis
   * - `error` is `std::cerr` with bold red emphasis
   */
  Printer_set();
};

//! \brief prints messages like "message... done"
class Task_message
{
  Printer& _printer;
  public:
  Task_message(Printer& p, std::string message, std::string sep = " ") : _printer(p) {_printer(message + "..." + sep);}
  ~Task_message() {_printer("done\n");}
};

}
#endif
