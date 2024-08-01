#include <hexed/Iges_parser.hpp>

namespace hexed {

Int Iges_parser::read_int(std::string) {
  return 0;
}

double Iges_parser::read_float(std::string) {
  return 0;
}

std::string Iges_parser::read_string(std::string) {
  return "";
}

Iges_parser::Iges_parser(std::string file_name) {
}

next::Sequence<const std::vector<std::string>&> Iges_parser::section(Section_id) {
  return {};
}

const std::vector<std::string>& Iges_parser::entry(Section_id, Int line) {
  return _entries[0][0];
}

}
