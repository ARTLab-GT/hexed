#ifndef HEXED_COMMAND_PARSER_HPP_
#define HEXED_COMMAND_PARSER_HPP_

#include <functional>
#include "Namespace.hpp"

namespace hexed
{

class Command_parser
{
  public:
  std::shared_ptr<Namespace> variables;
  std::map<std::string, std::function<void(std::string)>> statements;
  Command_parser();
  void exec(std::string commands);
};

}
#endif
