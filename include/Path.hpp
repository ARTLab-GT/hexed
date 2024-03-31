#ifndef HEXED_PATH_HPP_
#define HEXED_PATH_HPP_

#include <filesystem>
#include <vector>

namespace hexed
{

/*! \brief Finds a file in the working directory or system paths.
 * \details The "system paths" are:
 * - `/`
 * - `/local`
 * - `/usr`
 * - `usr/local`
 * - `$HOME`
 * - `$HOME/.local`
 */
class Path
{
  std::filesystem::path _home;
  std::vector<std::filesystem::path> _paths;

  public:
  using path = std::filesystem::path;
  //! \brief constructs a `Path` object with `subdir` appended to all search paths
  //! \details E.g., if the file your looking for is a shared library, you probably want to do `subdir = "lib"`
  Path(path subdir = {});
  //! \brief returns the list of system paths to search, including the `subdir` suffix
  inline std::vector<path> paths() {return _paths;}
  /*! \brief finds the file `target`
   * \details If `target` is a relative path, it will be searched for in the following directories are searched,
   * in order of precedence:
   * - All directories in `extra_dirs`, which defaults to the current working directory
   * - All directories in the environment variable `HEXEDPATH`, which should be a colon-separated list, if it is set
   * - All the system paths, defined above.
   *
   * If `target` is not found, the return value is an empty path.
   * If `target` is an absolute path, it will be returned if it exists and the empty path will be returned otherwise.
   */
  path find(path target, std::vector<path> extra_dirs = {{"."}});
};

}
#endif
