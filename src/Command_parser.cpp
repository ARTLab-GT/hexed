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
    if (std::isalpha(commands[place]) || commands[place] == '_') {
      std::string name = "";
      while (std::isalpha(commands[place]) || std::isdigit(commands[place])  || commands[place] == '_') {
        name.push_back(commands[place++]);
      }
      skip_spaces();
      HEXED_ASSERT(commands[place++] == '=', "expected assignment operator `=` after variable name");
      skip_spaces();
      std::string value;
      while (std::isdigit(commands[place])) value.push_back(commands[place++]);
      variables->assign(name, std::atoi(value.c_str()));
      HEXED_ASSERT(!more() || commands[place] == '\n', "expected end of line after assignment statement");
    }
    break;
  }
}

}
