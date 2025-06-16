#ifndef HEXED_PRINTER_HPP_
#define HEXED_PRINTER_HPP_

#include <iostream>
#include <vector>
#include <memory>

namespace hexed {

//! \brief abstract base class for handling different variations of printing things for user
//! \details basically a slightly higher-level version of output streams
class Printer {
  public:
  virtual ~Printer() = default;
  /*!
   * \param message what you want to print
   * \param emph If `true`, some form of emphasis formatting will be used, if possible.
   * \param replace If `true`, `message` will *replace* the current line of text on the screen,
   *                rather than be appended to it.
   */
  virtual void operator()(std::string message, bool emph = false, bool replace = false) = 0;
};

//! \brief `Printer` which is just a combination of other printers.
//! \details Anything you tell it to print will be forwarded to all of `Compound_printer::printers`
class Compound_printer : public Printer {
  public:
  std::vector<std::shared_ptr<Printer>> printers; //!< \brief list of printers---feel free to modify
  Compound_printer() = default;
  Compound_printer(std::vector<std::shared_ptr<Printer>>);
  void operator()(std::string message, bool emph = false, bool replace = false) override;
};

//! \brief prints to a `std::ostream`.
class Stream_printer : public Printer {
  std::ostream& _stream;
  std::string _format_code;
  std::string _reset_code;
  std::string _replace_code;

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
   * The latter will explicitly send an [ASCII escape code](https://en.wikipedia.org/wiki/ANSI_escape_code#Colors)
   * which sets the color to the default.
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
    //! \brief if `true`, any colors specified refer to the text background rather than the text itself
    bool background = false;
  };
  static const Format default_format;

  /*! \brief constructs a `Stream_printer`
   * \param stream stream to print to
   * \param use_escape_codes If `true`, [ASCII escape codes](https://en.wikipedia.org/wiki/ANSI_escape_code#Colors)
   *   will be used to apply formatting to the text and play other tricks.
   *   This is necessary for the `emph` and `replace` arguments of `operator()` to have any effect,
   *   but it also requires the stream to support ASCII escape codes.
   *   `cout` and `cerr` support these codes, but ASCII files do not.
   * \param emph_format Formatting to apply to emphasized text.
   *   Any value other than the default will require `use_escape_codes = true`.
   */
  Stream_printer(std::ostream& stream = std::cout, bool use_escape_codes = false, Format emph_format = default_format);

  void operator()(std::string, bool emph = false, bool replace = false) override;
};

/*! \brief Global `Printer` objects that any function can use to print messages.
 * \details These should be used instead of the builtin printing facilities,
 * so that they can be used to redirect the output of the entire program.
 * While making them global variables is conceptually ugly,
 * the convenience benefit over passing references to `Printer`s
 * to any code that might ever want to print something is substantial.
 */
namespace printers {
  extern Compound_printer info; //!< \brief general information
  extern Compound_printer warn; //!< \brief (non-fatal) warnings
  extern Compound_printer error; //!< \brief (fatal) errors
}

//! \brief prints messages like "message... done"
class Task_message {
  Printer& _printer;
  std::string _sep1;
  public:
  //! \details `sep0` inserted after ellipsis, `sep1` inserted before "done".
  Task_message(Printer& p, std::string message, std::string sep0 = " ", std::string sep1 = "");
  ~Task_message();
};

}
#endif
