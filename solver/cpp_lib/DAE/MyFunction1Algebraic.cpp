#include "MyFunction1Algebraic.hpp"
#include "Defs.hpp"
#include <math.h>
class MyFunction1Algebraic;

Real MyFunction1Algebraic::operator()(VariableArray& vaDifferential,
    VariableArray& vaAlgebraic, Real aCurrentTime) const
{
    return vaAlgebraic[1] * vaAlgebraic[1] 
        + vaDifferential[1] * vaDifferential[1] 
        - vaAlgebraic[2] - vaAlgebraic[0];
}

