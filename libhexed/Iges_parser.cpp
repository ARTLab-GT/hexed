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

Iges_parser::Iges_parser(std::string file_name) : _param_delim{0}, _record_delim{0}, _entries(5) {
  HEXED_ASSERT(std::filesystem::exists(file_name),
               format_str(1000, "`%s` is not an existing file", file_name.c_str()));
  std::ifstream file(file_name);
  std::vector<std::string> rec;
  std::string field;
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
        int i = 0;
        Section_id id;
        if (sec == 'G') {
          id = global;
          if (!_param_delim) {
            if (line[0] == ',') {
              _param_delim = ',';
              ++i;
            } else {
              HEXED_ASSERT(line[0] == '1' && line[1] == 'H', "failed to identify parameter delimiter");
              _param_delim = line[i + 2];
              i += 4;
            }
          }
          if (!_record_delim) {
            if (line[i] == _param_delim) {
              _record_delim = ';';
              ++i;
            } else {
              HEXED_ASSERT(line[i] == '1' && line[i + 1] == 'H', "failed to identify record delimiter");
              _record_delim = line[i + 2];
              HEXED_ASSERT(line[i + 3] == _param_delim, "no parameter delimiter after record delimiter specification");
              i += 4;
            }
          }
        } else {
          id = parameter;
          HEXED_ASSERT(_param_delim && _record_delim, "parameter section before delimiter specification");
        }
        int h_count = 0;
        for (; i < 72; ++i) {
          if ((line[i] == _param_delim || line[i] == _record_delim) && !h_count) {
            rec.push_back(field);
            field.clear();
            if (line[i] == _record_delim) {
              _entries[id].push_back(rec);
              rec.clear();
            }
          } else {
            h_count = std::max(0, h_count - 1);
            if (line[i] == 'H' && std::all_of(field.begin(), field.end(), [](char c){return std::isdigit(c);})) {
              h_count = std::stoi(field);
            }
            if (!field.empty() || line[i] != ' ') field.insert(field.size(), 1, line[i]);
          }
        }
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
