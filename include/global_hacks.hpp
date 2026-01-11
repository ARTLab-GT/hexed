#ifdef HEXED_USE_GLOBAL_HACKS
#ifndef HEXED_GLOBAL_HACKS_HPP_
#define HEXED_GLOBAL_HACKS_HPP_
#include <map>
#include <vector>
#include "Stopwatch_tree.hpp"

//! \namespace hexed::global_hacks
//! \brief Don't use any of this. This is not the file you're looking for. I'm sorry, but it's time for you to leave.
namespace hexed::global_hacks {

extern std::map<std::string, int> debug_message;
extern Stopwatch_tree stopwatch;
extern std::vector<double> numbers;

}
#endif
#endif
