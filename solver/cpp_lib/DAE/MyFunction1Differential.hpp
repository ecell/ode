#ifndef __MYFUNCTION1DIFFERENTIAL_HPP
#define __MYFUNCTION1DIFFERENTIAL_HPP

#include "Defs.hpp"
#include "Function.hpp"

class MyFunction1Differential
  :
  public Function
{
public:
    virtual ~MyFunction1Differential() {}

    virtual Real operator()(VariableArray& variable_array_differential,
        VariableArray& variable_array_algebraic, Real a_current_time) const;
};

#endif /* __MYFUNCTION1DIFFERENTIAL_HPP */
