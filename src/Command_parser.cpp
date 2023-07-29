#include <Command_parser.hpp>

namespace hexed
{

Command_parser::Command_parser() : variables{std::make_shared<Namespace>()} {}

void Command_parser::exec(std::string commands)
{
  int place = 0;
  auto more = [&](){return int(commands.size()) > place;};
  auto skip_spaces = [&](){while (commands[place] == ' ') ++place;}; // note: commands[commands.size()] exists and is a null character so this is safe
  while (more()) {
    skip_spaces();
    if (commands[place] == '\n') ++place;
    else if (std::isalpha(commands[place]) || commands[place] == '_') {
      std::string name = "";
      while (std::isalpha(commands[place]) || std::isdigit(commands[place])  || commands[place] == '_') {
        name.push_back(commands[place++]);
      }
      skip_spaces();
      HEXED_ASSERT(commands[place++] == '=', "expected assignment operator `=` after variable name");
      skip_spaces();
      HEXED_ASSERT(more(), "unexpected end of line in assignment statement");
      if (std::isdigit(commands[place]) || commands[place] == '.') {
        std::string value;
        bool is_int = true;
        while (std::isdigit(commands[place]) || commands[place] == '.' || std::tolower(commands[place]) == 'e') {
          is_int = is_int && std::isdigit(commands[place]);
          value.push_back(commands[place++]);
        }
        if (is_int) variables->assign(name, std::stoi(value.c_str()));
        else        variables->assign(name, std::stod(value.c_str()));
      } else {
        HEXED_ASSERT(false, format_str(100, "failed to parse value starting with `%c`", commands[place]));
      }
      skip_spaces();
      HEXED_ASSERT(!more() || commands[place] == '\n', "expected end of line after assignment statement");
    }
  }
}

}
