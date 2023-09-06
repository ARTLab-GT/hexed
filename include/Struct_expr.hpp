#ifndef HEXED_STRUCT_EXPR_HPP_
#define HEXED_STRUCT_EXPR_HPP_

#include "Interpreter.hpp"

namespace hexed
{

//! \brief evaluates a _structured expression_ in HIL
class Struct_expr
{
  public:
  const std::vector<std::string> names; //!< names of output variables
  const std::vector<std::string> exprs; //!< expressions for values of output variables
  /*! \anchor struct_expr
   * \param code A _structured expression_, meaning an expression that always assigns to the same numeric variables,
   *   regardless of the environment in which it was executed.
   *   Specifically, it must satisfy the following requirements
   *   - It consists only of assignment statements where the LHS is non-empty and does not contain the macro-substitution character `$`.
   *   - The RHS must not contain line breaks or semicolons (except the one that terminates the statement, of course).
   *   - The RHS of the assignment statements may contain `$`, but any resulting substitutions must not line breaks/semicolons or exit the interpreter.
   *   - Each RHS must evaluate to an `int` or `double`.
   */
  Struct_expr(std::string code);
  /*! \brief executes the code and returns the values of the output variables.
   * \details Variables are returned in the same order as `names`.
   * Evaluates the expressions with the provided interpreter.
   */
  std::vector<double> eval(Interpreter&) const;
};

}
#endif
