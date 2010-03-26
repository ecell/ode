#include "MyFunction1Differential.hpp"
#include "Defs.hpp"
#include <math.h>
class MyFunction1Differential;

Real MyFunction1Differential::operator()(VariableArray& vaDifferential,
    VariableArray& vaAlgebraic, Real aCurrentTime) const
{
    return -vaDifferential[0] * vaAlgebraic[2];
}

