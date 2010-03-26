#ifndef __MYFUNCTION2ALGEBRAIC_HPP
#define __MYFUNCTION2ALGEBRAIC_HPP

#include "Defs.hpp"
#include "Function.hpp"

class MyFunction2Algebraic
  :
  public Function
{
public:
    virtual ~MyFunction2Algebraic() {}

    virtual Real operator()(VariableArray& variable_array_differential, 
        VariableArray& variable_array_algebraic, Real a_current_time) const;
};

#endif /* __MYFUNCTION2AlGEBRAIC_HPP */
