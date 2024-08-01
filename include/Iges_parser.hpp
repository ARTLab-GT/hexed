#ifndef HEXED_IGES_PARSER_HPP_
#define HEXED_IGES_PARSER_HPP_

#include <string>
#include <fstream>
#include "math.hpp"
#include "Sequence.hpp"

namespace hexed {

class Iges_parser {
  public:
  enum Section_id {flag, start, global, directory, parameter, terminate};

  static Int read_int(std::string);
  static double read_float(std::string);
  static std::string read_string(std::string);

  Iges_parser(std::string file_name);
  next::Sequence<const std::vector<std::string>&> section(Section_id);
  const std::vector<std::string>& entry(Section_id, Int line);

  private:
  char _param_delim;
  char _record_delim;
  std::vector<std::vector<std::vector<std::string>>> _entries;
};

}
#endif
