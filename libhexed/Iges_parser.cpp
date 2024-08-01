#include <filesystem>
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

Iges_parser::Iges_parser(std::string file_name) : _entries(5) {
  HEXED_ASSERT(std::filesystem::exists(file_name),
               format_str(1000, "`%s` is not an existing file", file_name.c_str()));
  std::ifstream file(file_name);
  while (!file.eof()) {
    char line [81];
    file.getline(line, 81);
    int len = std::strlen(line);
    if (len != 0) {
      HEXED_ASSERT(len == 80, "line is shorter than 80 characters");
      char sec = line[72];
      HEXED_ASSERT(sec != 'C' && sec != 'B', "Only uncompressed ASCII IGES format is supported",
                   assert::Not_implemented_error)
      if (sec == 'S') {
        _entries[start].push_back({std::string(line, 72)});
      } else if (sec == 'G' || sec == 'P') {
      } else if (sec == 'D' || sec == 'T') {
      } else HEXED_THROW(format_str(100, "section character '%c' not recognized", sec));
    }
  }
}

next::Sequence<const std::vector<std::string>&> Iges_parser::section(Section_id sec) {
  return next::Sequence<const std::vector<std::string>&>::vector_view(_entries[sec]);
}

const std::vector<std::string>& Iges_parser::entry(Section_id sec, Int line) {
  HEXED_ASSERT(line > 0 && line <= Int(_entries[sec].size()), "line number out of bounds");
  return _entries[sec][line - 1];
}

}
