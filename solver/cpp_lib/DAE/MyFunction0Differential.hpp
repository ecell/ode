#ifndef __MYFUNCTION0DIFFERENTIAL_HPP
#define __MYFUNCTION0DIFFERENTIAL_HPP

#include "Defs.hpp"
#include "Function.hpp"

class MyFunction0Differential
  :
  public Function
{
public:
    virtual ~MyFunction0Differential() {}

    virtual Real operator()(VariableArray& variable_array_differential,
        VariableArray& variable_array_algebraic, Real a_current_time) const;
};

#endif /* __MYFUNCTION0DIFFERENTIAL_HPP */
