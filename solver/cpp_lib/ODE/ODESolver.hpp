#ifndef __ODESOLVER_HPP
#define __ODESOLVER_HPP

#ifdef DLL_EXPORT
#undef DLL_EXPORT
#include <gsl/gsl_rng.h>
#define DLL_EXPORT
#else
#include <gsl/gsl_rng.h>
#endif

#include <string>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "Defs.hpp"
#include "Function.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <boost/multi_array.hpp>

#include <gsl/gsl_poly.h>

typedef boost::multi_array<Real, 2> RealMatrix;

class ODESolver
{
    class Interpolant
    {      
    public:
        Interpolant(ODESolver *a_solver): 
            the_solver_(a_solver)
        {
            ; // do nothing
        }

        virtual ~Interpolant()
        {
        }
        virtual const Real get_difference(Time a_time, Real an_interval, 
            VariableArray::size_type the_index) const;

    private:
        ODESolver* the_solver_;
    };


public:
    ODESolver();
    virtual ~ODESolver();

    /** Register a function.
    * Functions must be registered in order of variables of array.
    * At least, one function must be set before step().
    */
    void register_function(Function* function)
    {
        the_function_.push_back(function);
    }

    /** Initialize the Solver
    * This method must be called before step().
    */
    void initialize(Real variable_array[], 
        const Integer a_variable_array_size);

    /** Get variables
    */
    void get_variable_array(Real variable_array[], 
        const Integer a_variable_array_size) const
    {
        std::copy(the_value_.begin(), the_value_.begin() + a_variable_array_size,
                  variable_array);
    }

    void initialize_solver(const Integer a_variable_array_size);

    /** Integrate variables.
    */
    virtual void integrate(Time a_time);

    /** Reschedule next step.
    * Return the next step time.
    */
    virtual Time reschedule();
    virtual void integrate_variable(Time a_time);

    /** Simulate one step.
    * Following three methods make a set.
    *  integrate()
    *  step()
    *  reschedule()
    */
    virtual Integer step();
    virtual Integer update_internal_state(Real a_step_interval);
    Real calculate_jacobian_norm();
    void calculate_jacobian();
    void initialize_radauIIA(const Integer the_system_size);
    void update_internal_state_radauIIA(Real a_step_interval);
    void update_internal_state_differential_stepper(Real a_step_interval);
    void set_variable_velocity(
        boost::detail::multi_array::sub_array<Real, 1> a_velocity_buffer);
    void set_jacobian_matrix(const Real a_step_interval);
    void decomp_jacobian_matrix();
    std::pair<bool, Real> calculate_radauIIA(const Real a_step_interval, 
        const Real a_previous_step_interval);
    void calculate_rhs(const Real a_step_interval);
    Real solve();
    Real estimate_local_error(const Real a_step_interval);
    virtual void update_internal_state_adaptive_differential_stepper(
        Real a_step_interval);
    void clear_variables();
    virtual bool calculate(Real a_step_interval);
    void inter_integrate();    
    void reset_all();
    virtual void reset();

    /** Register a status event.
    */
    void register_status_event(StatusEvent status_event)
    {
        the_status_event_.push_back(status_event);
    }
    Integer get_status_code();


    Real const get_number_equations()
    {
        return the_function_.size();
    }

    virtual const Integer get_order()
    {
        if (is_stiff_)
        {
            return 3;
        }
        else 
        {
            return 4;
        }
    }
  
    const Time get_current_time()
    {
        return the_current_time_;
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

            delete the_interpolant_;

            //error handling
            std::ostringstream ss1, ss2;

            ss1 << a_step_interval;
            std::string str1 = ss1.str();

            ss2 << the_min_step_interval_;
            std::string str2 = ss2.str();

            throw "The step interval (" + str1 
                + ") is behind the error-limit step interval (" + str2 + ").";
        }
        else
        {
            the_next_time_ = value;
        }
    }

    const Time get_next_time()
    {
        return the_next_time_;
    }


    Real get_tolerable_step_interval()
    {
        return the_tolerable_step_interval_;
    }
  
    const RealMatrix& get_taylor_series() const
    {
        return the_taylor_series_;
    }

    Real get_value(VariableArray::size_type an_index) const
    {
        return the_value_[an_index];
    }

private:
    Time the_current_time_;
    Time the_next_time_;
    Time the_last_time_;
  
    Real the_next_step_interval_;
    Real the_spectral_radius_;
    Real eta_, uround_;
    Real the_stopping_criterion_;
    Real the_tolerance_;
    Real rtoler_, atoler_;
    Real the_tolerable_step_interval_;
    Real alpha_, beta_, gamma_;
    Real the_jacobian_recalculate_theta_;
    Real the_absolute_tolerance_factor_;
    Real the_accepted_error_, the_accepted_step_interval_;
    Real the_max_error_ratio_;
    Real safety_;
    Real the_state_tolerance_factor_;
    Real the_derivative_tolerance_factor_;
    Real the_min_step_interval_;
    Real the_max_step_interval_;

    Integer check_interval_count_, switching_count_;
    Integer the_system_size_;
    Integer the_tolerable_rejected_step_count_;

    UnsignedInteger the_stiffness_counter_, the_rejected_step_counter_;
    UnsignedInteger the_max_iteration_number_;

    bool is_interrupted_;
    bool is_stiff_;  
    bool the_first_step_flag_, the_jacobian_calculate_flag_;
    bool the_state_flag_;

    RealVector the_eigen_vector_, the_temp_vector_;
    RealMatrix the_taylor_series_;
    RealMatrix the_jacobian_, the_w_;

    gsl_matrix* the_jacobian_matrix1_;
    gsl_permutation* the_permutation1_;
    gsl_vector* the_velocity_vector1_;
    gsl_vector* the_solution_vector1_;
  
    gsl_matrix_complex* the_jacobian_matrix2_;
    gsl_permutation* the_permutation2_;
    gsl_vector_complex* the_velocity_vector2_;
    gsl_vector_complex* the_solution_vector2_;
  
    VariableArray the_value_;
    VariableArray the_value_buffer_;
    FunctionArray the_function_;
    Interpolant* the_interpolant_;

    StatusEventArray the_status_event_;
};

#endif /* __ODESOLVER_HPP */
