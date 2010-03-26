#include <iostream>
#include <string>

#include "MyFunction0Differential.hpp"
#include "MyFunction1Differential.hpp"
#include "MyFunction0Algebraic.hpp"
#include "MyFunction1Algebraic.hpp"
#include "MyFunction2Algebraic.hpp"
#include "DAESolver.hpp"
#include "Defs.hpp"

void test_run(DAESolver& s, Time an_end_time)
{
    Time the_time(s.get_current_time());

    while (the_time < an_end_time)
    {
        s.integrate(the_time);

        //20100224追加
        Integer status_code = s.step();
        if (status_code != 0)
        {
            break;
        }

        std::cout <<"                           時刻とy:" << s.get_current_time()
            <<  " " << s.get_value_algebraic(0)<< std::endl;

        the_time = s.reschedule();
    }

}

int main()
{
    const Integer a_variable_array_size = 5;
    const Integer a_variable_array_size_algebraic = 3;
    const Integer a_variable_array_size_differential
        = a_variable_array_size - a_variable_array_size_algebraic;
    try
    {
        if (a_variable_array_size <= 0)
        {
            throw 1;
        }
        if (a_variable_array_size_algebraic <= 0)
        {
            throw 2;
        }
        if (a_variable_array_size_differential <= 0)
        {
            throw 3;
        }
    }
    catch (int e)
    {
        if (e == 1)
        {
            std::cout << "The number of variables (or equations) \
                'a_variable_array_size' must be positive." << std::endl;
        }
        else if (e == 2)
        {
            std::cout << "The number of Algebraic equations \
                'a_variable_array_size_algebraic' must be positive."
                << std::endl;
        }
        else if (e == 3)
        {
            std::cout << "The number of Diffrential equations \
                'a_variable_array_size_differential' must be positive."
                << std::endl;
        }
        exit(1);
    }

    DAESolver s;

    Real variable_array_differential[a_variable_array_size_differential]
        = {0.0, 1.0};
    Real variable_array_algebraic[a_variable_array_size_algebraic]
        = {-1.0, 0.0, 2.0};

    MyFunction0Differential f0_differential;
    MyFunction1Differential f1_differential;
    s.register_function(&f0_differential);
    s.register_function(&f1_differential);

    MyFunction0Algebraic f0_algebraic;
    MyFunction1Algebraic f1_algebraic;
    MyFunction2Algebraic f2_algebraic;
    s.register_function(&f0_algebraic);
    s.register_function(&f1_algebraic);
    s.register_function(&f2_algebraic);

    try
    {
        if (s.get_number_equations() != a_variable_array_size) throw 1;
    }
    catch (int e)
    {
        std::cout << "The number of equations must be equal to the \
            number of variables 'a_variable_array_size'."<< std::endl;
        exit(1);
    }

    const Time a_duration = 10.;
    try
    {
        if (a_duration <= 0.)
        {
            throw 1;
        }
    }
    catch (int e)
    {
        std::cout << "The time to end this simulation \
            'a_duration' must be positive." << std::endl;
        exit(1);
    }

    try
    {
        s.initialize(variable_array_differential, variable_array_algebraic,
            a_variable_array_size, a_variable_array_size_differential);

        //20100224追加
        StatusEvent se0 = {0, 0, 1.5, 1};
        StatusEvent se1 = {0, 1, 1.0, 2};
        StatusEvent se2 = {0, 1, 0.5, 3};
        StatusEvent se3 = {0, 0, 0.0, 4};
        StatusEvent se4 = {0, 1, -0.5, 5};
        StatusEvent se5 = {2, 1, 0.5, -6};

        s.register_status_event(se0);
        s.register_status_event(se1);
        s.register_status_event(se2);
        s.register_status_event(se3);
        s.register_status_event(se4);
        s.register_status_event(se5);

        test_run(s, a_duration);
    }
    catch (std::string str)
    {
        std::cout << str<<std::endl;
    }

    return 0;
}
