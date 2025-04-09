#include <hexed/utils.hpp>

namespace hexed {

std::string file_extension(std::string file_name) {
  unsigned extension_start = file_name.find_last_of(".");
  HEXED_ASSERT(extension_start != std::string::npos, "`file_name` has no extension");
  std::string case_sensitive = file_name.substr(extension_start + 1, std::string::npos);
  std::string ext = case_sensitive;
  for (char& c : ext) c = tolower(c);
  return ext;
}

std::string to_lower(std::string s) {
  for (char& c : s) c = std::tolower(c);
  return s;
}

}
