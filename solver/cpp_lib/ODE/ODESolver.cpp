#include <limits>
#include <cmath>

#include <boost/array.hpp>

#include <cstring>
#include <cstdlib>
#include <utility>
#include <cctype>
#include <functional>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>

#include "Defs.hpp"
#include "Function.hpp"
#include "ODESolver.hpp"

const Real SQRT6(2.4494897427831779);

ODESolver::ODESolver():
    the_current_time_(0.0),
    the_next_step_interval_(0.001),
    check_interval_count_(100),
    the_stiffness_counter_(0),
    is_stiff_(true),
    the_spectral_radius_(0.0),
    is_interrupted_(true),
    eta_(1.0),
    uround_(1e-10),
    the_stopping_criterion_(0.0),
    rtoler_(1e-6),
    atoler_(1e-6),
    the_tolerance_(1.0e-6),
    the_first_step_flag_(true),
    the_jacobian_calculate_flag_(true),
    the_system_size_(-1),
    the_jacobian_matrix1_(NULLPTR),
    the_permutation1_(NULLPTR),
    the_velocity_vector1_(NULLPTR),
    the_solution_vector1_(NULLPTR),
    the_jacobian_matrix2_(NULLPTR),
    the_permutation2_(NULLPTR),
    the_velocity_vector2_(NULLPTR),
    the_solution_vector2_(NULLPTR),
    switching_count_(20),
    the_tolerable_step_interval_(0.001),
    the_state_flag_(true),
    the_rejected_step_counter_(0),
    the_jacobian_recalculate_theta_(0.001),
    the_max_iteration_number_(7),
    the_absolute_tolerance_factor_(1.0),
    the_accepted_step_interval_(0.0),
    the_accepted_error_(0.0),
    the_tolerable_rejected_step_count_(std::numeric_limits< Integer >::max()),
    the_max_error_ratio_(1.0),
    the_state_tolerance_factor_(1.0),
    the_derivative_tolerance_factor_(1.0),
    safety_(0.9),
    the_next_time_(0.001),
    the_min_step_interval_(1e-100),
    the_max_step_interval_(std::numeric_limits<Real>::infinity()),
    the_last_time_(0.0)
{
    const Real pow913(pow(9.0, 1.0 / 3.0));
  
    alpha_ = (12.0 - pow913 * pow913 + pow913) / 60.0;
    beta_ = (pow913 * pow913 + pow913) * sqrt(3.0) / 60.0;
    gamma_ = (6.0 + pow913 * pow913 - pow913) / 30.0;
  
    const Real a_norm(alpha_ * alpha_ + beta_ * beta_);
  
    alpha_ /= a_norm;
    beta_ /= a_norm;
    gamma_ = 1.0 / gamma_;
  
    rtoler_ = 0.1 * pow(the_tolerance_, 2.0 / 3.0);
    atoler_ = rtoler_ * the_absolute_tolerance_factor_;

    the_value_.resize(0);
    the_value_buffer_.resize(0);
    the_eigen_vector_.resize(0);
    the_temp_vector_.resize(0);
    the_taylor_series_.resize(boost::extents[0]
        [static_cast<RealMatrix::index>(0)]);
    the_jacobian_.resize(boost::extents[0][0]);
    the_w_.resize(boost::extents[6][0]);
    the_function_.resize(0);
    the_status_event_.resize(0);

    the_interpolant_ = new Interpolant(this);
}

ODESolver::~ODESolver()
{
    if (the_jacobian_matrix1_)
    {
        gsl_matrix_free(the_jacobian_matrix1_);
    }
    if (the_permutation1_)
    {
        gsl_permutation_free(the_permutation1_);
    }
    if (the_velocity_vector1_)
    {
        gsl_vector_free(the_velocity_vector1_);
    }
    if (the_solution_vector1_)
    {
        gsl_vector_free(the_solution_vector1_);
    }
    if (the_jacobian_matrix2_)
    {
        gsl_matrix_complex_free(the_jacobian_matrix2_);
    }
    if (the_permutation2_)
    {
        gsl_permutation_free(the_permutation2_);
    }
    if (the_velocity_vector2_) 
    {
        gsl_vector_complex_free(the_velocity_vector2_);
    }
    if (the_solution_vector2_)
    {
        gsl_vector_complex_free(the_solution_vector2_);
    }

}

void ODESolver::initialize(Real variable_array[], 
    const Integer a_variable_array_size)
{  
    the_value_.resize(a_variable_array_size);
    memcpy(&the_value_[0], variable_array, 
        a_variable_array_size * sizeof(Real));

    the_value_buffer_ = the_value_;

    the_last_time_ = the_current_time_;


    is_interrupted_ = true;


    the_taylor_series_.resize(boost::extents[4]
        [static_cast<RealMatrix::index>(a_variable_array_size)]);

    //******************

    initialize_solver(a_variable_array_size);
}


void ODESolver::initialize_solver(const Integer a_variable_array_size)
{
    //warning !!!!!
    //以下のコードをそのままにすると
    //常にstiffになってしまう。
    //解決策：******************に持っていく。
    is_stiff_ = true;
    the_stiffness_counter_ = 0;
    //warning !!!!!

    const Integer a_size(a_variable_array_size);

    if (is_stiff_)
    {
        initialize_radauIIA(a_size);
    }

    if (a_variable_array_size != the_system_size_)
    {
        the_w_.resize(boost::extents[6][a_size]);
    }

    the_system_size_ = a_size;
}  

void ODESolver::initialize_radauIIA(Integer a_new_system_size)
{
    eta_ = 1.0;
    the_stopping_criterion_ =
        std::max(10.0 * uround_ / rtoler_, std::min(0.03, sqrt(rtoler_)));
  
    the_first_step_flag_ = true;
    the_jacobian_calculate_flag_ = true;

    if (a_new_system_size != the_system_size_)
    {
        the_jacobian_.resize(
            boost::extents[a_new_system_size][a_new_system_size]);
        the_eigen_vector_.resize(a_new_system_size);
        the_temp_vector_.resize(a_new_system_size);

        if (the_jacobian_matrix1_)
        {
            gsl_matrix_free(the_jacobian_matrix1_);
            the_jacobian_matrix1_ = 0;
        }

        if (a_new_system_size > 0)
        {
            the_jacobian_matrix1_ = 
                gsl_matrix_calloc(a_new_system_size, a_new_system_size);
        }

        if (the_permutation1_)
        {
            gsl_permutation_free(the_permutation1_);
            the_permutation1_ = 0;
        }
      
        if (a_new_system_size > 0)
        {
            the_permutation1_ = gsl_permutation_alloc(a_new_system_size);
        }

        if (the_velocity_vector1_)
        {
            gsl_vector_free(the_velocity_vector1_);
            the_velocity_vector1_ = 0;
        }
      
        if (a_new_system_size > 0)
        {
            the_velocity_vector1_ = gsl_vector_calloc(a_new_system_size);
        }

        if (the_solution_vector1_)
        {
            gsl_vector_free(the_solution_vector1_);
            the_solution_vector1_ = 0;
        }

        if (a_new_system_size > 0)
        {
            the_solution_vector1_ = gsl_vector_calloc(a_new_system_size);
        }

        if (the_jacobian_matrix2_)
        {
            gsl_matrix_complex_free(the_jacobian_matrix2_);
            the_jacobian_matrix2_ = 0;
        }
      
        if (a_new_system_size > 0)
        {
            the_jacobian_matrix2_ = gsl_matrix_complex_calloc(
                a_new_system_size, a_new_system_size);
        }

        if (the_permutation2_)
        {
            gsl_permutation_free(the_permutation2_);
            the_permutation2_ = 0;
        }
      
        if (a_new_system_size > 0)
        {
            the_permutation2_ = gsl_permutation_alloc(a_new_system_size);
        }

        if (the_velocity_vector2_)
        {
            gsl_vector_complex_free(the_velocity_vector2_);
            the_velocity_vector2_ = 0;
        }
      
        if (a_new_system_size > 0)
        {
            the_velocity_vector2_ = 
                gsl_vector_complex_calloc(a_new_system_size);
        }

        if (the_solution_vector2_)
        {
            gsl_vector_complex_free(the_solution_vector2_);
            the_solution_vector2_ = 0;
        }
      
        if (a_new_system_size > 0)
        {
            the_solution_vector2_ = 
                gsl_vector_complex_calloc(a_new_system_size);
        }
    }
}

  
void ODESolver::integrate(Time a_time)
{
    integrate_variable(a_time);
    Real const a_step_interval(the_next_time_ - the_current_time_);
    the_current_time_ = a_time;
    the_next_time_ = a_time + a_step_interval;
}

void ODESolver::integrate_variable(Time a_time)
{
    const Time a_current_time(a_time);
    const Real an_interval(a_current_time - the_last_time_);

    if (an_interval == 0.0)
    {
        return;
    }
  
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        const Real a_velocity_sum(
            the_interpolant_->get_difference(a_current_time, an_interval, c));

        the_value_[c] += a_velocity_sum;
    }

    the_last_time_ = a_current_time;
}


Time ODESolver::reschedule()
{
    const Time a_local_time(the_current_time_);
    const Time a_new_step_interval(get_step_interval());
  
    return a_new_step_interval + a_local_time;  
}



Integer ODESolver::step()
{
    return update_internal_state(the_next_step_interval_);
}



Integer ODESolver::update_internal_state(Real a_step_interval)
{
    // check the stiffness of this system by previous results
    if (check_interval_count_ > 0)
    {
        if (the_stiffness_counter_ % check_interval_count_ == 1)
        {
            if (is_stiff_)
            {
                the_spectral_radius_ = calculate_jacobian_norm();
            }

            const Real lambdah(the_spectral_radius_ * a_step_interval);
            if (is_stiff_ == (lambdah < 3.3 * 0.8))
            {
                if (the_stiffness_counter_ > 
                    check_interval_count_ * switching_count_)
                {
                    is_stiff_ = static_cast<bool>(!is_stiff_);
                    initialize_solver(the_system_size_);
                }
            }
            else
            {
                the_stiffness_counter_ = 1;
            }
        }
      
        ++the_stiffness_counter_;
    }
  


    if (is_stiff_)
    {
        update_internal_state_radauIIA(a_step_interval);
    }
    else
    {
        update_internal_state_adaptive_differential_stepper(a_step_interval);
    }

    return get_status_code();
}

Integer ODESolver::get_status_code()
{
    Integer status_code(0);
    Real minimum_occurrence_time(INF);
    const RealMatrix::size_type a_taylor_size(get_order());
    RealVector tmp_min;
    tmp_min.resize(0);

    for (StatusEventArray::size_type i(0); i < the_status_event_.size(); ++i)
    {
        Integer index = the_status_event_[i].variable_index;
        Real threshold = the_status_event_[i].threshold;
        Integer code = the_status_event_[i].status_code;
    
        Real a[a_taylor_size + 1];

        a[0] = the_value_[index] - threshold;
        a[1] = the_taylor_series_[0][index];
        Real power = 1.0;
        for (Integer j(0); j < a_taylor_size - 1; ++j)
        {
            power = power * the_tolerable_step_interval_;
            a[j + 2] = the_taylor_series_[j + 1][index] / power;
        }

        Real z[2 * a_taylor_size];
        gsl_poly_complex_workspace * w
            = gsl_poly_complex_workspace_alloc(a_taylor_size + 1);
        gsl_poly_complex_solve(a, a_taylor_size + 1, w, z);
        gsl_poly_complex_workspace_free(w);

        for (Integer k(0); k < a_taylor_size; k++)
        {       
            const Real solution(z[2 * k]);

            if (z[2 * k + 1] == 0.0 && solution >= 0. 
                && solution < the_next_time_ - the_current_time_)
            {
                
                if (solution < minimum_occurrence_time)
                {
                    minimum_occurrence_time = solution;
                    status_code = code;

                    tmp_min.push_back(minimum_occurrence_time);
                }
            }

        }
    }


    Integer count = 0;
    for (Integer n(0); n < tmp_min.size(); n++)
    {
        if (tmp_min[n] == minimum_occurrence_time)
        {
            count++;
        }
    }
    if (count > 1) throw "Some status events occurred at the same time.";
  
    return status_code;
}

void ODESolver::update_internal_state_adaptive_differential_stepper(
    Real a_step_interval)
{
    the_state_flag_ = false;

    clear_variables();

    Integer the_rejected_step_counter_(0);
    const Real max_error(the_max_error_ratio_);
  
    while (!calculate(a_step_interval))
    {
        if (++the_rejected_step_counter_ >= 
            the_tolerable_rejected_step_count_)
        {
            delete the_interpolant_;

            //error handling
            std::ostringstream ss;
            ss << the_tolerable_rejected_step_count_;
            std::string str = ss.str();
  
            throw "The times of rejections of step calculation \
                exceeded a maximum tolerable count (" + str + ").";
        }
      
        //shrink it if the error exceeds 110%
        a_step_interval = a_step_interval * safety_ 
            * std::pow(max_error, -1.0 / get_order());
    }
  
    the_tolerable_step_interval_ = a_step_interval;

    the_state_flag_ = true;

    // grow it if error is 50% less than desired
    Real a_new_step_interval(a_step_interval);
    if (max_error < 0.5)
    {
        a_new_step_interval = a_new_step_interval * safety_ 
            * std::pow(max_error, -1.0 / (get_order() + 1));
    }
  
    the_next_step_interval_ = 
        std::min(a_new_step_interval, the_max_step_interval_);
    update_internal_state_differential_stepper(a_step_interval);
}


bool ODESolver::calculate(Real a_step_interval)
{
    const Real eps_rel(the_tolerance_);
    const Real eps_abs(the_tolerance_ * the_absolute_tolerance_factor_);
    const Real a_y(the_state_tolerance_factor_);
    const Real a_dydt(the_derivative_tolerance_factor_);
    
    const Time a_current_time(the_current_time_);
    
    // ========= 1 ===========
  
    if (is_interrupted_)
    {
        inter_integrate();
        set_variable_velocity(the_taylor_series_[0]);
      
        for (VariableArray::size_type c(0); c < the_system_size_; ++c)
        {
            the_value_[c] = 
                the_taylor_series_[0][c] * (1.0 / 5.0) * a_step_interval 
                + the_value_buffer_[c];
        }
    }
    else
    {
        for (VariableArray::size_type c(0); c < the_system_size_; ++c)
        {
            // get k1
            the_taylor_series_[0][c] = the_w_[5][c];
  
            the_value_[c] =  
                the_taylor_series_[0][c] * (1.0 / 5.0) * a_step_interval
                + the_value_buffer_[c];
        }
    }

    // ========= 2 ===========
    the_current_time_ = a_current_time + a_step_interval * 0.2;
    inter_integrate();
    set_variable_velocity(the_w_[0]);
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        the_value_[c] = (the_taylor_series_[0][c] * (3.0 / 40.0) 
            + the_w_[0][c] * (9.0 / 40.0)) * a_step_interval 
            + the_value_buffer_[c];
    }
    
    // ========= 3 ===========
    the_current_time_ = a_current_time + a_step_interval * 0.3;
    inter_integrate();
    set_variable_velocity(the_w_[1]);

    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        the_value_[c] = (the_taylor_series_[0][c] * (44.0 / 45.0)
            - the_w_[0][c] * (56.0 / 15.0) 
            + the_w_[1][c] * (32.0 / 9.0)) * a_step_interval 
            + the_value_buffer_[c];
    }
    
    // ========= 4 ===========
    the_current_time_ = a_current_time + a_step_interval * 0.8;
    inter_integrate();
    set_variable_velocity(the_w_[2]);
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        the_value_[c] = (the_taylor_series_[0][c] * (19372.0 / 6561.0)
            - the_w_[0][c] * (25360.0 / 2187.0)
            + the_w_[1][c] * (64448.0 / 6561.0)
            - the_w_[2][c] * (212.0 / 729.0)) * a_step_interval
            + the_value_buffer_[c];
    }

    // ========= 5 ===========
    the_current_time_ = a_current_time + a_step_interval * (8.0 / 9.0);
    inter_integrate();
    set_variable_velocity(the_w_[3]);
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        // temporarily set Y^6
        the_taylor_series_[1][c] = 
            the_taylor_series_[0][c] * (9017.0 / 3168.0)
            - the_w_[0][c] * (355.0 / 33.0)
            + the_w_[1][c] * (46732.0 / 5247.0)
            + the_w_[2][c] * (49.0 / 176.0)
            - the_w_[3][c] * (5103.0 / 18656.0);

        the_value_[c] = the_taylor_series_[1][c] * a_step_interval
            + the_value_buffer_[c];
    }
    
    // ========= 6 ===========
    
    // estimate stiffness
    Real a_denominator(0.0);
    Real a_spectral_radius(0.0);
    
    the_current_time_ = a_current_time + a_step_interval;
    inter_integrate();
    set_variable_velocity(the_w_[4]);

    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        the_taylor_series_[2][c] = the_taylor_series_[0][c] * (35.0 / 384.0)
            // + the_w_[0][c] * 0.0
            + the_w_[1][c] * (500.0 / 1113.0)
            + the_w_[2][c] * (125.0 / 192.0)
            + the_w_[3][c] * (-2187.0 / 6784.0)
            + the_w_[4][c] * (11.0 / 84.0);

        a_denominator +=
            (the_taylor_series_[2][c] - the_taylor_series_[1][c])
            * (the_taylor_series_[2][c] - the_taylor_series_[1][c]);

        the_value_[c] = 
            the_taylor_series_[2][c] * a_step_interval + the_value_buffer_[c];
    }
    
    // ========= 7 ===========
    the_current_time_ = a_current_time + a_step_interval;
    inter_integrate();
    set_variable_velocity(the_w_[5]);

    // evaluate error
    Real max_error(0.0);
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        // calculate error
        const Real an_estimated_error(
            (the_taylor_series_[0][c] * (71.0 / 57600.0)
            + the_w_[1][c] * (-71.0 / 16695.0)
            + the_w_[2][c] * (71.0 / 1920.0)
            + the_w_[3][c] * (-17253.0 / 339200.0)
            + the_w_[4][c] * (22.0 / 525.0)
            + the_w_[5][c] * (-1.0 / 40.0)) * a_step_interval);
      
        a_spectral_radius +=
            (the_w_[5][c] - the_w_[4][c]) * (the_w_[5][c] - the_w_[4][c]);
      
        // calculate velocity for Xn+.5
        the_taylor_series_[1][c] =
            the_taylor_series_[0][c] * (6025192743.0 / 30085553152.0)
            + the_w_[1][c] * (51252292925.0 / 65400821598.0)
            + the_w_[2][c] * (-2691868925.0 / 45128329728.0)
            + the_w_[3][c] * (187940372067.0 / 1594534317056.0)
            + the_w_[4][c] * (-1776094331.0 / 19743644256.0)
            + the_w_[5][c] * (11237099.0 / 235043384.0);
      
        const Real a_tolerance(
            eps_rel * (a_y * fabs(the_value_buffer_[c]) 
            + a_dydt * fabs(the_taylor_series_[2][c]) * a_step_interval) 
            + eps_abs);
        const Real an_error(fabs(an_estimated_error / a_tolerance));

        if (an_error > max_error)
        {
            max_error = an_error;
        }
    }

    a_spectral_radius /= a_denominator;
    a_spectral_radius = sqrt(a_spectral_radius);

    reset_all(); // reset all value

    the_max_error_ratio_ = max_error;
    the_current_time_ = a_current_time;
    
    if (max_error > 1.1)
    {
        // reset the stepper current time
        reset();
        is_interrupted_ = true;
        return false;
    }
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        const Real k1(the_taylor_series_[0][c]);
        const Real v_2(the_taylor_series_[1][c]);
        const Real v1(the_taylor_series_[2][c]);
        const Real k7(the_w_[5][c]);
      
        the_taylor_series_[1][c] = 
            -4.0 * k1 + 8.0 * v_2 - 5.0 * v1 + k7;
        the_taylor_series_[2][c] = 
            5.0 * k1 - 16.0 * v_2 + 14.0 * v1 - 3.0 * k7;
        the_taylor_series_[3][c] = 
            -2.0 * k1 + 8.0 * v_2 - 8.0 * v1 + 2.0 * k7;
    }
    
    // set the error limit interval
    is_interrupted_ = false;
    
    the_spectral_radius_ = a_spectral_radius / a_step_interval;

    return true;
}



void ODESolver::reset_all()
{
  the_value_ = the_value_buffer_;
}

void ODESolver::reset()
{
    // is this needed?
    for (RealMatrix::index i(0); i < 4; ++i)
    {
        for (RealMatrix::index j(0); j < the_system_size_; ++j)
        {
            the_taylor_series_[i][j] = 0.0;
        }
    }

    the_value_ = the_value_buffer_;
}



void ODESolver::inter_integrate()
{
    const Time a_current_time(the_current_time_);
    const Real an_interval(a_current_time - the_last_time_);

    if (an_interval > 0.0)
    {
        for (VariableArray::size_type c(0); c < the_system_size_; ++c)
        {
            const Real a_velocity_sum(
                the_interpolant_->get_difference(
                    a_current_time, an_interval, c));
            the_value_[c] += a_velocity_sum; 
        }
    }

}




void ODESolver::update_internal_state_radauIIA(Real a_step_interval)
{
    if (!the_jacobian_matrix1_)
    {
        update_internal_state_differential_stepper(a_step_interval);
        return;
    }

    the_state_flag_ = false;

    const Real a_previous_step_interval(get_step_interval());

    clear_variables();

    the_rejected_step_counter_ = 0;

    set_variable_velocity(the_w_[3]);

    if (the_jacobian_calculate_flag_ || is_interrupted_)
    {
        calculate_jacobian();
        set_jacobian_matrix(a_step_interval);
    }
    else
    {
        if (a_previous_step_interval != a_step_interval)
        {
            set_jacobian_matrix(a_step_interval);
        }
    }

    for (;;)
    {
        std::pair<bool, Real> const aResult(
            calculate_radauIIA(a_step_interval, a_previous_step_interval));

        if (aResult.first)
        {
            break;
        }
        a_step_interval = aResult.second;

        if (++the_rejected_step_counter_ >= 
            the_tolerable_rejected_step_count_)
        {
            delete the_interpolant_;

            std::ostringstream ss;
            ss << the_tolerable_rejected_step_count_;
            std::string str = ss.str();
    
            throw "The times of rejections of step calculation \
                exceeded a maximum tolerable count (" + str + ").";
        }

        if (!the_jacobian_calculate_flag_)
        {
            calculate_jacobian();
            the_jacobian_calculate_flag_ = true;
        }

        set_jacobian_matrix(a_step_interval);
    }
    
    the_tolerable_step_interval_ = a_step_interval;
    
    the_spectral_radius_ = calculate_jacobian_norm();
    
    // the_w_ will already be transformed to Z-form
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        the_w_[3][c] = the_w_[2][c];
        the_w_[3][c] /= a_step_interval;
        the_value_[c] = the_value_buffer_[c];
    }
    
    for (VariableArray::size_type c(0); c < the_system_size_; c++)
    {
        const Real z1(the_w_[0][c]);
        const Real z2(the_w_[1][c]);
        const Real z3(the_w_[2][c]);

        the_taylor_series_[0][c] = (13.0 + 7.0 * SQRT6) / 3.0 * z1
            + (13.0 - 7.0 * SQRT6) / 3.0 * z2
            + 1.0 / 3.0 * z3;
        the_taylor_series_[1][c] = -(23.0 + 22.0 * SQRT6) / 3.0 * z1
            + (-23.0 + 22.0 * SQRT6) / 3.0 * z2
            - 8.0 / 3.0 * z3;
        the_taylor_series_[2][c] = (10.0 + 15.0 * SQRT6) / 3.0 * z1
            + (10.0 - 15.0 * SQRT6) / 3.0 * z2
            + 10.0 / 3.0 * z3;

        the_taylor_series_[0][c] /= a_step_interval;
        the_taylor_series_[1][c] /= a_step_interval;
        the_taylor_series_[2][c] /= a_step_interval;
    }

    the_state_flag_ = true;
    
    // an extra calculation for resetting the activities of processes
    update_internal_state_differential_stepper(a_step_interval);
}

std::pair<bool, Real> ODESolver::calculate_radauIIA(
    const Real a_step_interval, const Real a_previous_step_interval)
{
    Real a_new_step_interval;
    Real a_norm;
    Real theta(fabs(the_jacobian_recalculate_theta_));
    
    UnsignedInteger an_iterator(0);
    
    if (!is_interrupted_)
    {
        const Real c1((4.0 - SQRT6) / 10.0);
        const Real c2((4.0 + SQRT6) / 10.0);
        const Real c3q(a_step_interval / a_previous_step_interval);
        const Real c1q(c3q * c1);
        const Real c2q(c3q * c2);
        for (VariableArray::size_type c(0); c < the_system_size_; ++c)
        {
            const Real cont3(the_taylor_series_[2][c]);
            const Real cont2(the_taylor_series_[1][c] + 3.0 * cont3);
            const Real cont1(the_taylor_series_[0][c] + 2.0 * cont2 
                - 3.0 * cont3);
    
            const Real z1(a_previous_step_interval * c1q *
                (cont1 + c1q * (cont2 + c1q * cont3)));
            const Real z2(a_previous_step_interval * c2q * 
                (cont1 + c2q * (cont2 + c2q * cont3)));
            const Real z3(a_previous_step_interval * c3q * 
                (cont1 + c3q * (cont2 + c3q * cont3)));
  
            the_w_[0][c] = 4.3255798900631553510 * z1
                + 0.33919925181580986954 * z2 + 0.54177053993587487119 * z3;
            the_w_[1][c] = -4.1787185915519047273 * z1
                - 0.32768282076106238708 * z2 + 0.47662355450055045196 * z3;
            the_w_[2][c] = -0.50287263494578687595 * z1
                + 2.5719269498556054292 * z2 - 0.59603920482822492497 * z3;
        }
    }
    else
    {
        for (VariableArray::size_type c(0); c < the_system_size_; ++c)
        {
            the_w_[0][c] = 0.0;
            the_w_[1][c] = 0.0;
            the_w_[2][c] = 0.0;
        }
    }
  
    eta_ = pow(std::max(eta_, uround_), 0.8);
  
    for (;;)
    {
        if (an_iterator == the_max_iteration_number_)
        {
            // XXX: this will be addressed somehow.
            // std::cerr << "matrix is repeatedly singular" << std::endl;
            break;
        }
      
        calculate_rhs(a_step_interval);
        const Real previous_norm(std::max(a_norm, uround_));
        a_norm = solve();

        if (an_iterator > 0 && (an_iterator != the_max_iteration_number_ - 1))
        {
            const Real a_theta_q = a_norm / previous_norm;
            if (an_iterator > 1)
            {
                theta = sqrt(a_theta_q * theta);
            }
            else
            {
                theta = a_theta_q;
            }

            if (theta < 0.99)
            {
                eta_ = theta / (1.0 - theta);
                const Real an_iteration_error(eta_ * a_norm * pow(theta, 
                    static_cast<int>(
                        the_max_iteration_number_ - 2 - an_iterator)) 
                    / the_stopping_criterion_);
                if (an_iteration_error >= 1.0)
                {
                    Real base = std::max(1e-4, 
                        std::min(20.0, an_iteration_error));
                    Real exponent = -1.0 / 
                        (4 + the_max_iteration_number_ - 2 - an_iterator);
                    return std::make_pair(false, a_step_interval * 0.8 * 
                        pow(base, exponent));
                }
            }
            else
            {
                return std::make_pair(false, a_step_interval * 0.5);
            }
        }

        if (eta_ * a_norm <= the_stopping_criterion_)
        {
            break;
        }
      
        an_iterator++;
    }
  
    // the_w_ is transformed to Z-form
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        const Real w1(the_w_[0][c]);
        const Real w2(the_w_[1][c]);
        const Real w3(the_w_[2][c]);
      
        the_w_[0][c] = w1 * 0.091232394870892942792
            - w2 * 0.14125529502095420843
            - w3 * 0.030029194105147424492;
        the_w_[1][c] = w1 * 0.24171793270710701896
            + w2 * 0.20412935229379993199
            + w3 * 0.38294211275726193779;
        the_w_[2][c] = w1 * 0.96604818261509293619 + w2;
    }
    
    const Real an_error(estimate_local_error(a_step_interval));
    
    Real a_safety_factor(
        std::min(0.9, 0.9 * (1 + 2 * the_max_iteration_number_) / 
        (an_iterator + 1 + 2 * the_max_iteration_number_)));
    a_safety_factor = 
        std::max(0.125, std::min(5.0, pow(an_error, 0.25) / a_safety_factor));
  
    a_new_step_interval = a_step_interval / a_safety_factor;
    

    if (an_error < 1.0)
    {
        // step is accepted
      
        if (!the_first_step_flag_)
        {
            Real a_gustafsson_factor(
                the_accepted_step_interval_ / a_step_interval * 
                pow(an_error * an_error / the_accepted_error_, 0.25) / 0.9);
            a_gustafsson_factor = 
                std::max(0.125, std::min(5.0, a_gustafsson_factor));
    
            if (a_safety_factor < a_gustafsson_factor)
            {
                a_safety_factor = a_gustafsson_factor;
                a_new_step_interval = a_step_interval / a_gustafsson_factor;
            }
        }

        the_accepted_step_interval_ = a_step_interval;
        the_accepted_error_ = std::max(1.0e-2, an_error);

        if (the_rejected_step_counter_ != 0)
        {
            a_new_step_interval = 
                std::min(a_new_step_interval, a_step_interval);
        }

        the_first_step_flag_ = false;

        const Real a_step_interval_rate(
            a_new_step_interval / a_step_interval);

        if (theta <= the_jacobian_recalculate_theta_)
        {
            the_jacobian_calculate_flag_ = false;
        }
        else
        {
            the_jacobian_calculate_flag_ = true;
        }

        if (a_step_interval_rate >= 1.0 && a_step_interval_rate <= 1.2)
        {
            the_next_step_interval_ = 
                std::min(a_step_interval, the_max_step_interval_);
        }
        else
        {
            the_next_step_interval_ = 
                std::min(a_new_step_interval, the_max_step_interval_);
        }

        return std::make_pair(true, a_step_interval);
    }
    else
    {
        // step is rejected
      
        if (the_first_step_flag_)
        {
            a_new_step_interval = 0.1 * a_step_interval;
        }
      
        return std::make_pair(false, a_new_step_interval);
    }
}

Real ODESolver::estimate_local_error(const Real a_step_interval)
{
    Real an_error;
  
    const Real hee1((-13.0 - 7.0 * SQRT6) / (3.0 * a_step_interval));
    const Real hee2((-13.0 + 7.0 * SQRT6) / (3.0 * a_step_interval));
    const Real hee3(-1.0 / (3.0 * a_step_interval));
    
    // the_w_ will already be transformed to Z-form
    for (VariableArray::size_type c(0); c < the_system_size_; c++)
    {
        gsl_vector_set(the_velocity_vector1_, c,
            the_w_[3][c] 
            + the_w_[0][c] * hee1
            + the_w_[1][c] * hee2
            + the_w_[2][c] * hee3);
    }
    
    gsl_linalg_LU_solve(the_jacobian_matrix1_, the_permutation1_,
        the_velocity_vector1_, the_solution_vector1_);
    
    an_error = 0.0;
    Real a_difference;
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        const Real a_tolerance(
            rtoler_ * fabs(the_value_buffer_[c]) + atoler_);
        a_difference = gsl_vector_get(the_solution_vector1_, c);

        // for the case (an_error >= 1.0)
        the_value_[c] = the_value_buffer_[c] + a_difference;

        a_difference /= a_tolerance;
        an_error += a_difference * a_difference;
    }
    
    an_error = std::max(sqrt(an_error / the_system_size_), 1e-10);
    
    if (an_error < 1.0)
    {
        return an_error;
    }
    
    if (the_first_step_flag_ || the_rejected_step_counter_ != 0)
    {
        set_variable_velocity(the_w_[4]);

        for (VariableArray::size_type c(0); c < the_system_size_; ++c)
        {
            gsl_vector_set(the_velocity_vector1_, c,
                the_w_[4][c]
                + the_w_[0][c] * hee1
                + the_w_[1][c] * hee2
                + the_w_[2][c] * hee3);
        }
      
        gsl_linalg_LU_solve(the_jacobian_matrix1_, the_permutation1_,
            the_velocity_vector1_, the_solution_vector1_);

        an_error = 0.0;
        for (VariableArray::size_type c(0); c < the_system_size_; ++c)
        {
            const Real a_tolerance(
                rtoler_ * fabs(the_value_buffer_[c]) + atoler_);
    
            Real a_difference(gsl_vector_get(the_solution_vector1_, c));
    
            a_difference /= a_tolerance;
  
            an_error += a_difference * a_difference;
        }

        an_error = std::max(sqrt(an_error / the_system_size_), 1e-10);
    }
  
    return an_error;
}


Real ODESolver::solve()
{
    gsl_linalg_LU_solve(the_jacobian_matrix1_, the_permutation1_,
        the_velocity_vector1_, the_solution_vector1_);
    gsl_linalg_complex_LU_solve(the_jacobian_matrix2_, the_permutation2_,
        the_velocity_vector2_, the_solution_vector2_);
  
    Real a_norm(0.0);
    Real delta_w(0.0);
    gsl_complex comp;
  
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        Real a_tolerance2(rtoler_ * fabs(the_value_buffer_[c]) + atoler_);
        a_tolerance2 = a_tolerance2 * a_tolerance2;
      
        delta_w = gsl_vector_get(the_solution_vector1_, c);
        the_w_[0][c] += delta_w;
        a_norm += delta_w * delta_w / a_tolerance2;
      
        comp = gsl_vector_complex_get(the_solution_vector2_, c);
      
        delta_w = GSL_REAL(comp);
        the_w_[1][c] += delta_w;
        a_norm += delta_w * delta_w / a_tolerance2;
      
        delta_w = GSL_IMAG(comp);
        the_w_[2][c] += delta_w;
        a_norm += delta_w * delta_w / a_tolerance2;
    }
  
    return sqrt(a_norm / (3 * the_system_size_));
}




void ODESolver::calculate_rhs(const Real a_step_interval)
{
    const Real alphah(alpha_ / a_step_interval);
    const Real betah(beta_ / a_step_interval);
    const Real gammah(gamma_ / a_step_interval);
  
    const Time a_current_time(the_current_time_);

    gsl_complex comp;
    
    RealVector tif (the_system_size_ * 3);
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        const Real z(the_w_[0][c] * 0.091232394870892942792
            - the_w_[1][c] * 0.14125529502095420843
            - the_w_[2][c] * 0.030029194105147424492);
      
        the_value_[c] = the_value_buffer_[c] + z;
    }
  
    // ========= 1 ===========
  
    the_current_time_ = 
        a_current_time + a_step_interval * (4.0 - SQRT6) / 10.0;
    set_variable_velocity(the_w_[4]);
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        tif[c] = the_w_[4][c] * 4.3255798900631553510;
        tif[c + the_system_size_] = the_w_[4][c] * -4.1787185915519047273;
        tif[c + the_system_size_ * 2] = 
            the_w_[4][c] * -0.50287263494578687595;

        const Real z(the_w_[0][c] * 0.24171793270710701896
            + the_w_[1][c] * 0.20412935229379993199
            + the_w_[2][c] * 0.38294211275726193779);

        the_value_[c] = the_value_buffer_[c] + z;
    }
  
    // ========= 2 ===========
  
    the_current_time_ = a_current_time 
        + a_step_interval * (4.0 + SQRT6) / 10.0;
    set_variable_velocity(the_w_[4]);
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        tif[c] += the_w_[4][c] * 0.33919925181580986954;
        tif[c + the_system_size_] -= the_w_[4][c] * 0.32768282076106238708;
        tif[c + the_system_size_ * 2] += the_w_[4][c] * 2.5719269498556054292;
      
        const Real z(the_w_[0][c] * 0.96604818261509293619 + the_w_[1][c]);

        the_value_[c] = the_value_buffer_[c] + z;
    }
    
    // ========= 3 ===========
  
    the_current_time_ = a_current_time + a_step_interval;
    set_variable_velocity(the_w_[4]);
    
    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        tif[c] += the_w_[4][c] * 0.54177053993587487119;
        tif[c + the_system_size_] += the_w_[4][c] * 0.47662355450055045196;
        tif[c + the_system_size_ * 2] -= 
            the_w_[4][c] * 0.59603920482822492497;
      
        const Real w1(the_w_[0][c]);
        const Real w2(the_w_[1][c]);
        const Real w3(the_w_[2][c]);
      
        gsl_vector_set(the_velocity_vector1_, c, tif[c] - w1 * gammah);

        const Real p1(tif[c + the_system_size_] - w2 * alphah + w3 * betah);
        const Real p2(tif[c + the_system_size_ * 2] 
            - w2 * betah - w3 * alphah);
        GSL_SET_COMPLEX(&comp, p1, p2);

        gsl_vector_complex_set(the_velocity_vector2_, c, comp);
    }
    

    the_current_time_ = a_current_time;
}





void ODESolver::set_jacobian_matrix(const Real a_step_interval)
{
    const Real alphah(alpha_ / a_step_interval);
    const Real betah(beta_ / a_step_interval);
    const Real gammah(gamma_ / a_step_interval);
  
    gsl_complex comp1, comp2;
  
    for (RealVector::size_type i(0); i < the_system_size_; i++)
    {
        for (RealVector::size_type j(0); j < the_system_size_; j++)
        {
            const Real a_partial_derivative(the_jacobian_[i][j]);
            gsl_matrix_set(the_jacobian_matrix1_, i, j, a_partial_derivative);
  
            GSL_SET_COMPLEX(&comp1, a_partial_derivative, 0);
            gsl_matrix_complex_set(the_jacobian_matrix2_, i, j, comp1);
        }
    }
  
    for (VariableArray::size_type c(0); c < the_system_size_; c++)
    {
        const Real a_partial_derivative(
            gsl_matrix_get(the_jacobian_matrix1_, c, c));
        gsl_matrix_set(the_jacobian_matrix1_, c, c, 
            gammah + a_partial_derivative);
      
        comp1 = gsl_matrix_complex_get(the_jacobian_matrix2_, c, c);
        GSL_SET_COMPLEX(&comp2, alphah, betah);
        gsl_matrix_complex_set(the_jacobian_matrix2_, c, c,
            gsl_complex_add(comp1, comp2));
    }
  
    decomp_jacobian_matrix();
}


void ODESolver::decomp_jacobian_matrix()
{
    int a_sign_num;
  
    gsl_linalg_LU_decomp(the_jacobian_matrix1_, the_permutation1_, 
        &a_sign_num);
    gsl_linalg_complex_LU_decomp(the_jacobian_matrix2_, the_permutation2_, 
        &a_sign_num);
}

void ODESolver::update_internal_state_differential_stepper(
    Real a_step_interval)
{  
    set_step_interval(a_step_interval);

    // check if the step interval was changed, by epsilon
    if (std::fabs(the_tolerable_step_interval_ - a_step_interval)
       > std::numeric_limits<Real>::epsilon())
    {
        is_interrupted_ = true;
    }  
    else
    {
        is_interrupted_ = false;
    }
}


void ODESolver::set_variable_velocity(
    boost::detail::multi_array::sub_array<Real, 1> a_velocity_buffer)
{  
    for (RealMatrix::index i(0);
        i < static_cast<RealMatrix::index>(a_velocity_buffer.size()); ++i)
    {
        a_velocity_buffer[i] = 0.0;
    }  

    for (FunctionArray::size_type i(0); i < the_system_size_; i++)
    {
        a_velocity_buffer[i] = 
            (*the_function_[i])(the_value_, the_current_time_);
    } 
}



void ODESolver::calculate_jacobian()
{
    Real a_perturbation;
  
    for (VariableArray::size_type i(0); i < the_system_size_; ++i)
    {
        const Real a_value(the_value_[i]);
      
        a_perturbation = sqrt(uround_ * std::max(1e-5, fabs(a_value)));
        the_value_[i] = the_value_buffer_[i] + a_perturbation;
      
        set_variable_velocity(the_w_[4]);

        for (VariableArray::size_type j(0); j < the_system_size_; ++j)
        {
            the_jacobian_[j][i] = 
                -(the_w_[4][j] - the_w_[3][j]) / a_perturbation;
        }
      
        the_value_[i] = a_value;
    }

}



Real ODESolver::calculate_jacobian_norm()
{
    std::fill(the_eigen_vector_.begin(), the_eigen_vector_.end(), 
        sqrt(1.0 / the_system_size_));

    Real sum, norm;

    for (Integer count(0); count < 3; count++)
    {
        norm = 0.0;
        for (RealVector::size_type i(0); i < the_system_size_; i++)
        {
            sum = 0.0;
            for (RealVector::size_type j(0); j < the_system_size_; j++)
            {
                const Real a_partial_derivative(the_jacobian_[i][j]);
                sum += a_partial_derivative * the_eigen_vector_[j];
            }
            the_temp_vector_[i] = sum;
  
            norm += the_temp_vector_[i] * the_temp_vector_[i];
        }

        norm = sqrt(norm);
      
        for (RealVector::size_type i(0); i < the_system_size_; i++)
        {
            the_eigen_vector_[i] = the_temp_vector_[i] / norm;
        }
    }
  
    return norm;
}

void ODESolver::clear_variables()
{
  // save original value values
  the_value_buffer_ = the_value_;
}



const Real ODESolver::Interpolant::get_difference(
    Time a_time, Real an_interval, VariableArray::size_type the_index) const
{
    if (!the_solver_->the_state_flag_)
    {
        return 0.0;
    }

    const Time a_current_time(the_solver_->get_current_time());
    const Real a_time_interval1(a_time - a_current_time);
    const Real a_time_interval2(a_time_interval1 - an_interval);
  
    const RealMatrix& a_taylor_series(the_solver_->get_taylor_series());
    const Real* a_taylor_coefficient_ptr(a_taylor_series.origin() + the_index);
  

    // calculate first order.
    // here it assumes that always a_taylor_size >= 1
    Real a_value1(*a_taylor_coefficient_ptr * a_time_interval1);
    Real a_value2(*a_taylor_coefficient_ptr * a_time_interval2);

    // check if second and higher order calculations are necessary.
    const RealMatrix::size_type a_taylor_size(the_solver_->get_order());
    if (a_taylor_size >= 2)
    {
        const Real a_step_interval_inv(
            1.0 / the_solver_->get_tolerable_step_interval());
      
        const RealMatrix::size_type a_stride(a_taylor_series.strides()[0]);
      
        Real a_factorial_inv1(a_time_interval1);
        Real a_factorial_inv2(a_time_interval2);
      
        RealMatrix::size_type s(a_taylor_size - 1);
      
        const Real theta1(a_time_interval1 * a_step_interval_inv);
        const Real theta2(a_time_interval2 * a_step_interval_inv);
      
        do 
        {
            // main calculation for the 2+ order
            a_taylor_coefficient_ptr += a_stride;
            const Real a_taylor_coefficient(*a_taylor_coefficient_ptr);
          
            a_factorial_inv1 *= theta1;
            a_factorial_inv2 *= theta2;
          
            a_value1 += a_taylor_coefficient * a_factorial_inv1;
            a_value2 += a_taylor_coefficient * a_factorial_inv2;
          
            --s;
        } while (s != 0);
    }
  
    return a_value1 - a_value2;
}


