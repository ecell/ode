#ifndef __MYFUNCTION1ALGEBRAIC_HPP
#define __MYFUNCTION1ALGEBRAIC_HPP

#include "Defs.hpp"
#include "Function.hpp"

class MyFunction1Algebraic
  :
  public Function
{
public:
    virtual ~MyFunction1Algebraic() {}

    virtual Real operator()(VariableArray& variable_array_differential,
        VariableArray& variable_array_algebraic, Real a_current_time) const;
};

#endif /* __MYFUNCTION1AlGEBRAIC_HPP */
