#ifndef __DAESOLVER_HPP
#define __DAESOLVER_HPP

#include <string>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include "Defs.hpp"
#include "Function.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <boost/multi_array.hpp>

//20100224追加
#include <gsl/gsl_poly.h>

const Real SQRT6(2.4494897427831779);

typedef boost::multi_array<Real, 2> RealMatrix;

class DAESolver
{
    class Interpolant
    {
    public:
        Interpolant(DAESolver *a_solver)
          :
          the_solver_(a_solver)
        {
            ; // do nothing
        }

        virtual ~Interpolant()
        {
        }
        virtual const Real get_difference(Time a_time,
                    Real an_interval,
                    VariableArray::size_type the_index) const;

    private:
        DAESolver* the_solver_;
    };

public:
    DAESolver();
    virtual ~DAESolver();

    /** Integrate variables.
     */
    virtual void integrate(Time a_time);

    virtual void integrate_variable(Time a_time);

    /** Register a function.
     *
     * Functions must be registered in order of variables of array.
     * Functions of differential must be set first, 
     * then functions of algebraic must be set.
     * At least, one function must be set before step().
     */
    void register_function(Function* fun)
    {
        the_function_.push_back(fun);
    }

    /** Initialize the solver
     *
     * This method must be called before step().
     */
    void initialize(Real variable_array_differential[],
        Real variable_array_algebraic[], 
        const Integer a_variable_array_size,
        const Integer a_variable_array_size_differential);

    /** Get variables of differential
     */
    void get_variable_differential_array(Real variables[],
        const Integer a_variable_array_size) const
    {
        memcpy(variables, &the_value_differential_[0],
            a_variable_array_size * sizeof(Real));
    }

    /** Get variables of algebraic
     */
    void get_variable_algebraic_array(Real variables[],
        const Integer a_variable_array_size) const
    {
        memcpy(variables, &the_value_algebraic_[0],
            a_variable_array_size * sizeof(Real));
    }

    virtual Integer update_internal_state(Real a_step_interval);
    void clear_variables();

    /** Simulate one step.
     *
     * Following three methods make a set.
     *  integrate()
     *  step()
     *  reschedule()
     */
    virtual Integer step();

    void set_variable_velocity(
        boost::detail::multi_array::sub_array<Real, 1> a_velocity_buffer);
    void calculate_jacobian();
    void set_jacobian_matrix(Real const a_step_interval);
    void decomp_jacobian_matrix();
    std::pair< bool, Real > calculate_radauIIA(
        Real const a_step_interval, Real const a_previous_step_interval);
    void calculate_rhs(Real a_step_interval);
    Real solve();
    Real estimate_local_error(Real const a_step_interval);
    virtual void update_internal_state_differential_stepper(
        Real a_step_interval);

    /** Resuchedule next step.
     *
     * Reteurn the next step time.
     */
    virtual Time reschedule();

    //20100224追加
    /** Register a status event.
     */
    void register_status_event(StatusEvent se)
    {
        the_status_event_.push_back(se);
    }

    Integer get_status_code();

    Real const get_number_equations()
    {
        return the_function_.size();
    }

    const Time get_current_time()
    {
        return the_current_time_;
    }
    const RealMatrix & get_taylor_series() const
    {
        return the_taylor_series_;
    }
    Real get_tolerable_step_interval()
    {
        return the_tolerable_step_interval_;
    }

    const Real get_step_interval()
    {
        const Time a_next_time(the_next_time_);
        const Time a_current_time(the_current_time_);
        if (a_current_time == INF)
        {
            return INF;
        }

        return a_next_time - a_current_time;
    }

    void set_step_interval(const Real value)
    {
        set_next_time(the_current_time_ + value);
    }

    void set_next_time(Real value)
    {
        const Real a_step_interval(value - the_current_time_);

        if (a_step_interval < the_min_step_interval_)
        {
            the_next_time_ = value;

            //error handling
            std::ostringstream ss1, ss2;

            ss1 << a_step_interval;
            std::string str1 = ss1.str();

            ss2 << the_min_step_interval_;
            std::string str2 = ss2.str();
	
            delete the_interpolant_;

            throw "The step interval ("
                + str1 
                + ") is behind the error-limit step interval ("
                + str2 + ").";
        }
        else
        {
            the_next_time_ = value;
        }
    }

    Real get_value_differential(VariableArray::size_type the_index) const
    {
        return the_value_differential_[the_index];
    }

    Real get_value_algebraic(VariableArray::size_type the_index) const
    {
        return the_value_algebraic_[the_index];
    }

private:
    Time the_current_time_;
    Time the_next_time_;
    Time the_last_time_;

    Real the_tolerable_step_interval_;
    Real eta_, uround_;
    Real the_stopping_criterion_;
    Real alpha_, beta_, gamma_;
    Real the_absolute_tolerance_, atoler_;
    Real the_relative_tolerance_, rtoler_;
    Real the_next_step_interval_;
    Real the_jacobian_recalculate_theta_;
    Real the_accepted_error_, the_accepted_step_interval_;
    Real the_min_step_interval_;
    Real the_max_step_interval_;

    Integer the_system_size_;
    Integer the_function_differential_size_;

    UnsignedInteger the_max_iteration_number_;

    bool the_state_flag_;
    bool the_first_step_flag_;
    bool the_jacobian_calculate_flag_;
    bool is_interrupted_;
    bool the_rejected_step_flag_;

    RealMatrix the_taylor_series_;

    VariableArray the_value_differential_;
    VariableArray the_value_algebraic_;
    VariableArray the_value_differential_buffer_;
    VariableArray the_value_algebraic_buffer_;

    FunctionArray the_function_;

    gsl_matrix* the_jacobian_matrix1_;
    gsl_permutation* the_permutation1_;
    gsl_vector* the_velocity_vector1_;
    gsl_vector* the_solution_vector1_;

    gsl_matrix_complex* the_jacobian_matrix2_;
    gsl_permutation* the_permutation2_;
    gsl_vector_complex* the_velocity_vector2_;
    gsl_vector_complex* the_solution_vector2_;

    Interpolant* the_interpolant_;

    std::vector<RealVector> the_jacobian_;

    RealVector the_w_;
    RealVector the_activity_algebraic_buffer_;

    //20100224追加
    StatusEventArray the_status_event_;
};

#endif /* __DAESOLVER_HPP */
