#include <iostream>
#include "cpp_lib/ODE/ODESolver.hpp"
#include "cpp_lib/DAE/DAESolver.hpp"

typedef struct _StepResult
{
    long int state_event;
    int err;
}StepResult;

StepResult call_step(ODESolver *solver)
{
    StepResult step_result = {0, 0};

    try
    {
        step_result.state_event = solver->step();
    }
    catch (char* error_message)
    {
        std::cout << error_message << std::endl;
        step_result.err = 1;
    }

    return step_result;
}

StepResult call_step(DAESolver *solver)
{
    StepResult step_result = {0, 0};

    try
    {
        step_result.state_event = solver->step();
    }
    catch (char* error_message)
    {
        std::cout << error_message << std::endl;
        step_result.err = 1;
    }

    return step_result;
}

