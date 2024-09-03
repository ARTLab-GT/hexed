#ifndef HEXED_IGES_PARSER_HPP_
#define HEXED_IGES_PARSER_HPP_

#include <string>
#include <fstream>
#include "math.hpp"
#include "Sequence.hpp"

namespace hexed {

/*! \brief Parses IGES files to obtain parameter values.
 * \details Parses files in the
 * [Initial Graphics Exchange Specification, Version 6](https://filemonger.com/specs/igs/devdept.com/version6.pdf)
 * format.
 * This class does not interpret the geometric information in the files;
 * that is done by the `brep::Geom` class.
 * It simply reads the file and breaks it down into a list of _sections_ (e.g., _Start_, _Global_, _Directory_, etc.),
 * where each section contains a list of _entries_, and each entry is a list of strings representing data _fields_.
 * The fields can then be interpreted as integers, floating-point values,
 * or [Hollerith strings](https://en.wikipedia.org/wiki/Hollerith_constant)
 * with `read_int()`, `read_float()`, and `read_string()`.
 */
class Iges_parser {
  public:
  //! \brief Identifies the sections of an IGES file
  enum Section_id {start, global, directory, parameter, terminate};

  //! \brief Reads an integer value from a human-readable string
  //! \details An empty string defaults to 0.
  static Int read_int(std::string);

  //! \brief Reads a floating-point value from a human-readable string
  static double read_float(std::string);

  //! \brief Translates a Hollerith string to a normal string.
  //! \details E.g., `"6HRuffin"` translates to `"Ruffin"`
  static std::string read_string(std::string);

  /*! \brief Constructs a parser and reads a file.
   * \details `file_name` must be a path, absolute or relative, to an IGES file
   * with an \ref add_geom "appropriate extension".
   */
  Iges_parser(std::string file_name);

  //! \brief get the list of entries for a given section
  const next::Sequence<const std::vector<std::string>&> section(Section_id) const;

  /*! \brief gets an entry by which section it is in and its line number in the file
   * \details In IGES, a "pointer" to an entity
   * is an integer indicating the line number of its entry in the _Directory_ section.
   * Thus, this function can be used to resolve pointers.
   */
  const std::vector<std::string>& entry(Section_id, Int line) const;

  private:
  char _param_delim;
  char _record_delim;
  std::vector<std::vector<std::vector<std::string>>> _entries;
  std::vector<std::vector<Int>> _line_map;
};

}
#endif
