#include "MyFunction0Differential.hpp"
#include "Defs.hpp"
#include <math.h>
class MyFunction0Differential;

Real MyFunction0Differential::operator()(
    VariableArray& variable_array_differential,
    VariableArray& variable_array_algebraic, Real a_current_time) const
{
    return variable_array_differential[1];
}

