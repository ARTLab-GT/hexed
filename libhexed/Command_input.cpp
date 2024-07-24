#include <iostream>
#include <cstdio>
#include <termios.h>
#include <unistd.h>
#include <Command_input.hpp>
#include <assert.hpp>
#include <utils.hpp>

namespace hexed
{

const int Command_input::unlimited = -1;

Command_input::Command_input(int n_history) : _n_hist{n_history} {}

void read_char(char* c, int n)
{
  HEXED_ASSERT(read(fileno(stdin), c, n) == n, "failed to read keyboard input");
}

std::string Command_input::get()
{
  // add entry to history buffer
  while (_history.size() > size_t(_n_hist)) _history.pop_back();
  _history.emplace_front();
  auto display = _history.begin();
  auto begin = _history.begin();
  // set up termios
  termios old_settings;
  tcgetattr(fileno(stdin), &old_settings);
  termios new_settings = old_settings;
  new_settings.c_lflag &= ~ECHO;
  new_settings.c_lflag &= ~ICANON;
  tcsetattr(fileno(stdin), TCSANOW, &new_settings);
  // get current cursor position
  std::cout << "\x1b[6n" << std::flush;
  char response [100] {};
  for (char* data = response;; ++data) {
    read_char(data, 1);
    if (*data == 'R') break;
  }
  HEXED_ASSERT(response[0] == 27 && response[1] == '[', "attempt to obtain cursor position received unintelligible response");
  int line_start = std::stoi(response + std::string(response).find(';') + 1);
  // function to set the current command to a new (modified) command rather than an unmodified entry from the history
  auto modify = [&]() {
    if (display != begin) {
      *begin = *display;
      display = begin;
    }
  };
  // obtain input
  // There will be a string (the latest history entry) and a cursor position.
  // Every time the user enters a character, first the input string and the cursor position variables are updated,
  // and then the display is updated to match their current values.
  char c;
  int pos = 0; // cursor position
  do {
    read_char(&c, 1); // get character
    if (std::isprint(c)) { // if printable, add to input line
      modify();
      begin->insert(pos++, 1, c);
    } else if (c == 127) { // if backspace, delete the previous character
      if (pos) {
        --pos;
        modify();
        begin->erase(pos, 1);
      }
    } else if (c == 27) { // if escape sequence...
      std::string escape = "  ";
      read_char(escape.data(), 2); // get rest of escape sequence
      if (escape == "[A") { // if up arrow, switch to previous history entry
        if (display < _history.end() - 1) {
          ++display;
          pos = display->size();
        }
      } else if (escape == "[B") { // if down arrow, switch to next history entry
        if (display > begin) {
          --display;
          pos = display->size();
        }
      } else if (escape == "[C") { // if right arrow, move cursor right
        if (pos < int(display->size())) ++pos;
      } else if (escape == "[D") { // if left arrow, move cursor left
        if (pos) --pos;
      } else HEXED_ASSERT(false, format_str(100, "could not parse escape sequence `%s` in keyboard input", escape.c_str()));
    }
    // update the text displayed on the screen
    printf("\x1b[%iG\x1b[K%s\x1b[%iG", line_start, display->c_str(), line_start + pos);
    std::cout << std::flush;
  } while (c != '\n' && c != EOF);
  std::cout << std::endl;
  // delete latest history entry if it doesn't contain new information
  std::string input = *display;
  modify();
  if (begin->size() == 0) {
    _history.pop_front();
  } else if (_history.size() > 1) if (_history[0] == _history[1]) _history.pop_front();
  // clean up termios
  tcsetattr(fileno(stdin), TCSANOW, &old_settings);
  return input;
}

}
