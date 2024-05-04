#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <Command_input.hpp>
#include <assert.hpp>
#include <utils.hpp>

namespace hexed
{

int Command_input::_instances = 0;
termios Command_input::_old_settings;
const int Command_input::unlimited = -1;

Command_input::Command_input(int n_history)
: _n_hist{n_history}
{
  if (!_instances++) {
    tcgetattr(fileno(stdin), &_old_settings);
    termios new_settings = _old_settings;
    new_settings.c_lflag &= ~ECHO;
    new_settings.c_lflag &= ~ICANON;
    tcsetattr(fileno(stdin), TCSANOW, &new_settings);
  }
}

Command_input::Command_input(const Command_input& other)
: _n_hist{other._n_hist}, _history{other._history}
{
  ++_instances;
}

Command_input::~Command_input() {if (!--_instances) tcsetattr(fileno(stdin), TCSANOW, &_old_settings);}

void read_char(char* c, int n)
{
  HEXED_ASSERT(read(fileno(stdin), c, n) == n, "failed to read keyboard input");
}

std::string Command_input::get()
{
  _history.emplace_front();
  auto display = _history.begin();
  auto begin = _history.begin();
  char c;
  int pos = 0;
  auto modify = [&]() {
    if (display != begin) {
      *begin = *display;
      display = begin;
    }
  };
  do {
    read_char(&c, 1);
    if (std::isprint(c)) {
      modify();
      begin->insert(pos++, 1, c);
    } else if (c == 127) {
      if (pos) {
        --pos;
        modify();
        begin->erase(pos, 1);
      }
    } else if (c == 27) {
      std::string escape = "  ";
      read_char(escape.data(), 2);
      if (escape == "[A") {
        if (display < _history.end() - 1) {
          ++display;
          pos = display->size();
        }
      } else if (escape == "[B") {
        if (display > begin) {
          --display;
          pos = display->size();
        }
      } else if (escape == "[C") {
        if (pos < int(display->size())) ++pos;
      } else if (escape == "[D") {
        if (pos) --pos;
      } else HEXED_ASSERT(false, format_str(100, "could not parse escape sequence `%s` in keyboard input", escape.c_str()));
    }
    std::cout << "\x1b[2K\x1b[0G" << *display << std::string(display->size() - pos, '\b') << std::flush;;
  } while (c != '\n' && c != EOF);
  std::cout << std::endl;
  std::string input = *display;
  modify();
  if (begin->size() == 0) {
    _history.pop_front();
  } else if (_history.size() > 1) if (_history[0] == _history[1]) _history.pop_front();
  return input;
}

}
