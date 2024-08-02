#include <filesystem>
#include <map>
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

Iges_parser::Iges_parser(std::string file_name) : _param_delim{0}, _record_delim{0}, _entries(5), _line_map(5) {
  HEXED_ASSERT(std::filesystem::exists(file_name),
               format_str(1000, "`%s` is not an existing file", file_name.c_str()));
  std::ifstream file(file_name);
  std::vector<std::string> rec;
  std::string field;
  std::map<char, Section_id> _section_chars {
    {'S', start},
    {'G', global},
    {'D', directory},
    {'P', parameter},
    {'T', terminate},
  };
  while (!file.eof()) {
    char line [81];
    file.getline(line, 81);
    int len = std::strlen(line);
    if (len != 0) {
      HEXED_ASSERT(len == 80, "line is shorter than 80 characters");
      HEXED_ASSERT(line[72] != 'C' && line[72] != 'B', "Only uncompressed ASCII IGES format is supported",
                   assert::Not_implemented_error)
      HEXED_ASSERT(_section_chars.count(line[72]), format_str(100, "section character '%c' not recognized", line[72]));
      Section_id sec = _section_chars[line[72]];
      _line_map[sec].push_back(_entries[sec].size());
      if (sec == start) {
        _entries[start].push_back({std::string(line, 72)});
      } else if (sec == global || sec == parameter) {
        int i = 0;
        if (sec == global) {
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
        } else HEXED_ASSERT(_param_delim && _record_delim, "parameter section before delimiter specification");
        int h_count = 0;
        for (; i < 64 + 8*(sec == global); ++i) {
          if ((line[i] == _param_delim || line[i] == _record_delim) && !h_count) {
            rec.push_back(field);
            field.clear();
            if (line[i] == _record_delim) {
              _entries[sec].push_back(rec);
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
      } else if (sec == directory || sec == terminate) {
        for (int i_block = 0; i_block < 9; ++i_block) {
          int start = 8*i_block;
          while (line[start] == ' ' && start < 8*(i_block + 1)) ++start;
          rec.emplace_back(line + start, line + 8*(i_block + 1));
        }
        if (sec == terminate || rec.size() == 18) {
          _entries[sec].push_back(rec);
          rec.clear();
        }
      }
    }
  }
}

next::Sequence<const std::vector<std::string>&> Iges_parser::section(Section_id sec) {
  return next::Sequence<const std::vector<std::string>&>::vector_view(_entries[sec]);
}

const std::vector<std::string>& Iges_parser::entry(Section_id sec, Int line) {
  HEXED_ASSERT(line > 0 && line <= Int(_line_map[sec].size()), "line number out of bounds");
  int entry = _line_map[sec][line - 1];
  HEXED_ASSERT(entry < Int(_entries.size()), "error mapping lines to entries");
  return _entries[sec][entry];
}

}
