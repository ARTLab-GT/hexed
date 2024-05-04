#ifndef HEXED_COMMAND_INPUT_HPP_
#define HEXED_COMMAND_INPUT_HPP_

#include <termios.h>
#include <string>
#include <deque>

namespace hexed
{

class Command_input
{
  static int _instances;
  static termios _old_settings;
  int _n_hist;
  std::deque<std::string> _history;
  public:
  static const int unlimited;
  Command_input(int n_history = unlimited);
  Command_input(const Command_input&);
  ~Command_input();
  std::string get();
};

}

#endif
