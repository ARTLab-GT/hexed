#ifndef HEXED_STRUCT_EXPR_HPP_
#define HEXED_STRUCT_EXPR_HPP_

#include "Interpreter.hpp"

namespace hexed
{

//! \brief evaluates a _structured expression_ in HIL
class Struct_expr
{
  Interpreter& _inter;
  public:
  const std::vector<std::string> names; //!< names of output variables
  const std::vector<std::string> exprs; //!< expressions for values of output variables
  /*!
   * \param inter `code` will be evaluated in a \ref Interpreter::make_sub "sub-intepreter" of `inter`.
   * \param code A _structured expression_, meaning an expression that always assigns to the same numeric variables,
   *   regardless of the environment in which it was executed.
   *   Specifically, it must satisfy the following requirements
   *   - It consists only of assignment statements where the LHS is non-empty and does not contain the macro-substitution character `$`.
   *   - Multiline string literals are forbidden (assign them to a global variable or use `newline` from the standard library).
   *   - The RHS of the assignment statements may contain `$`, but any resulting substitutions must not add statements or exit the interpreter.
   *   - Each RHS must evaluate to an `int` or `double`.
   */
  Struct_expr(Interpreter& inter, std::string code);
  std::vector<double> eval(); //!< \brief executes the code and returns the values of the output variables in the same order as `names`
};

}
#endif
