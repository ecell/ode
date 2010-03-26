#include "MyFunction0.hpp"
#include "Defs.hpp"
#include <math.h>
class MyFunction0;

Real MyFunction0::operator()(VariableArray& va, Real a_current_time) const
{
    return -va[0];
}

