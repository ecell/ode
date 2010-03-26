#ifndef __MYFUNCTION0ALGEBRAIC_HPP
#define __MYFUNCTION0ALGEBRAIC_HPP

#include "Defs.hpp"
#include "Function.hpp"

class MyFunction0Algebraic
  :
  public Function
{
public:
    virtual ~MyFunction0Algebraic() {}

    virtual Real operator()(VariableArray& variable_array_differential,
        VariableArray& variable_array_algebraic, Real a_current_time) const;
};

#endif /* __MYFUNCTION0AlGEBRAIC_HPP */
