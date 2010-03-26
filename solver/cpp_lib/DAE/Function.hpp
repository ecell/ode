
#ifndef __FUNCTION_HPP
#define __FUNCTION_HPP

#include "Defs.hpp"

class Function
{
public:
    virtual ~Function() {}

    virtual Real operator()(VariableArray& variable_array_differential,
        VariableArray& variable_array_algebraic, 
        Real a_current_time) const = 0;
};

DECLARE_TYPE(std::vector<Function*>, FunctionArray);

#endif /* __FUNCTION_HPP */
