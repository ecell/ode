#ifndef __MYFUNCTION0_HPP
#define __MYFUNCTION0_HPP

#include "Defs.hpp"
#include "Function.hpp"

class MyFunction0 
  :
  public Function
{
public:
    virtual ~MyFunction0() {}
  
    virtual Real operator()(VariableArray& va, Real a_current_time) const;
};

#endif /* __MYFUNCTION0_HPP */
