#include "MyFunction2Algebraic.hpp"
#include "Defs.hpp"
#include <math.h>
class MyFunction2Algebraic;

Real MyFunction2Algebraic::operator()(
    VariableArray& variable_array_differential,
    VariableArray& variable_array_algebraic, Real a_current_time) const
{
    return variable_array_differential[0] * variable_array_differential[0]
        + variable_array_algebraic[0] * variable_array_algebraic[0] - 1.0;
}

