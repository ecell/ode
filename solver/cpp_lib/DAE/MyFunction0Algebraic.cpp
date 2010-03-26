#include "MyFunction0Algebraic.hpp"
#include "Defs.hpp"
#include <math.h>
class MyFunction0Algebraic;

Real MyFunction0Algebraic::operator()(
    VariableArray& variable_array_differential,
    VariableArray& variable_array_algebraic, Real a_current_time) const
{
    return 2.0 * variable_array_differential[0]
        * variable_array_differential[1]
        + 2.0 * variable_array_algebraic[0] * variable_array_algebraic[1];
}

