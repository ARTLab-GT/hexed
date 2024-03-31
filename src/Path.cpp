#include <cstdlib>
#include <config.hpp>
#include <Path.hpp>

namespace hexed
{

Path::Path(path subdir)
: _home{std::getenv("HOME")}, _paths{"/", "/local", "/usr", "/usr/local", _home, _home/".local"}
{
  for (path& p : _paths) p /= subdir;
}

Path::path Path::find(path target, std::vector<path> extra_dirs)
{
  // add pathes in `$HEXEDPATH`
  const char* env_hexedpath = std::getenv("HEXEDPATH");
  if (env_hexedpath) {
    std::string hexedpath = env_hexedpath;
    std::string::size_type start = 0;
    while (start < hexedpath.size()) {
      std::string::size_type end = std::min(hexedpath.find(":", start), hexedpath.size());
      extra_dirs.emplace_back(hexedpath.substr(start, end - start));
      start = end + 1;
    }
  }
  // search for `target`.
  // Note that search proceeds in reverse order, but if `target` is found later, it will override previous value.
  path found;
  for (std::vector<path>* dirs : {&_paths, &extra_dirs}) {
    for (path p : *dirs) {
      path candidate = p/target;
      if (std::filesystem::exists(candidate)) {
        found = candidate;
      }
    }
  }
  return found;
}

}
