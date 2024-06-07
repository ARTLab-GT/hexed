#ifndef HEXED_COMMAND_INPUT_HPP_
#define HEXED_COMMAND_INPUT_HPP_

#include <string>
#include <deque>

namespace hexed
{

//! \brief Handles command line input with history.
class Command_input
{
  int _n_hist;
  std::deque<std::string> _history;
  public:
  static const int unlimited; //!< \brief opaque value used to communicate no limit on the size of the history buffer
  /*! \brief Constructs a `Command_input` and sets the history buffer size.
   * \details All calls to `get()` with this `Command_input` instance will share the same history.
   * If `n_history` is not `Command_input::unlimited`,
   * whenever the history reaches a size of `n_history + 1` entries,
   * the earliest entry will be forgotten to keep a size of `n_history`.
   */
  Command_input(int n_history = unlimited);
  /*! \brief Obtains a line of input from the command line.
   * \details Waits for input until a newline or EOF is received.
   * This input is then returned as a string and added to the history buffer.
   * The terminal newline or EOF is not included, and only printable characters are included.
   * The up and down arrows can be used to access and navigate the history buffer.
   * Empty lines are not added to the history, nor are any lines that are identical to the previous line.
   */
  std::string get();
};

}

#endif
