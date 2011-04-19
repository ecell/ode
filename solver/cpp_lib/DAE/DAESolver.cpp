#include <limits>
#include <cmath>

#include <boost/array.hpp>

#include <vector>
#include <string>
#include <sstream>

#include "Defs.hpp"
#include "Function.hpp"
#include "DAESolver.hpp"

#include <iostream>
#include <assert.h>
#include <iomanip>

DAESolver::DAESolver()
  :
  the_next_time_(0.001),
  the_current_time_(0.0),
  the_last_time_(0.0),
  the_system_size_(-1),
  the_function_differential_size_(-1),
  the_jacobian_matrix1_(NULLPTR),
  the_permutation1_(NULLPTR),
  the_velocity_vector1_(NULLPTR),
  the_solution_vector1_(NULLPTR),
  the_jacobian_matrix2_(NULLPTR),
  the_permutation2_(NULLPTR),
  the_velocity_vector2_(NULLPTR),
  the_solution_vector2_(NULLPTR),
  the_state_flag_(true),
  the_tolerable_step_interval_(0.001),
  eta_(1.0),
  the_stopping_criterion_(0.0),
  uround_(1e-10),
  the_absolute_tolerance_(1e-6),
  the_relative_tolerance_(1e-6),
  the_first_step_flag_(true),
  the_jacobian_calculate_flag_(true),
  is_interrupted_(true),
  the_next_step_interval_(0.001),
  the_rejected_step_flag_(false),
  the_jacobian_recalculate_theta_(0.001),
  the_max_iteration_number_(7),
  the_accepted_step_interval_(0.0),
  the_accepted_error_(0.0),
  the_min_step_interval_(1e-100),
  the_max_step_interval_(std::numeric_limits<Real>::infinity())
{
    const Real pow913(pow(9.0, 1.0 / 3.0));
  
    alpha_ = (12.0 - pow913 * pow913 + pow913) / 60.0;
    beta_ = (pow913 * pow913 + pow913) * sqrt(3.0) / 60.0;
    gamma_ = (6.0 + pow913 * pow913 - pow913) / 30.0;
    
    const Real a_norm(alpha_ * alpha_ + beta_ * beta_);
    
    alpha_ /= a_norm;
    beta_ /= a_norm;
    gamma_ = 1.0 / gamma_;
    
    const Real a_ratio(the_absolute_tolerance_ / the_relative_tolerance_);
    rtoler_ = 0.1 * pow(the_relative_tolerance_, 2.0 / 3.0);
    atoler_ = rtoler_ * a_ratio;

    the_interpolant_ = new Interpolant(this);

    the_value_differential_.resize(0);
    the_value_algebraic_.resize(0);
    the_value_differential_buffer_.resize(0);
    the_value_algebraic_buffer_.resize(0);

    the_function_.resize(0);

    the_taylor_series_.resize(boost::extents[0]
        [static_cast< RealMatrix::index >(0)]);

    the_jacobian_.resize(0);
    the_w_.resize(0);
    the_activity_algebraic_buffer_.resize(0);
    the_status_event_.resize(0);
}

DAESolver::~DAESolver()
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

void DAESolver::integrate(Time a_time)
{
    integrate_variable(a_time);
    Real const a_step_interval(the_next_time_ - the_current_time_);
    the_current_time_ = a_time;
    the_next_time_ = a_time + a_step_interval;
}

void DAESolver::integrate_variable(Time a_time)
{
    const Time a_current_time(a_time);
    const Real an_interval(a_current_time - the_last_time_);

    if (an_interval == 0.0)
    {
        return;
    }
  
    for (VariableArray::size_type c(0);
        c < the_function_differential_size_; ++c)
    {
        const Real a_velocity_sum(
            the_interpolant_->get_difference(a_current_time, an_interval, c));

        the_value_differential_[c] += a_velocity_sum;
    }

    for (VariableArray::size_type c(the_function_differential_size_);
        c < the_system_size_; c++)
    {
        const Real a_velocity_sum(
            the_interpolant_->get_difference(a_current_time, an_interval, c));
        const VariableArray::size_type an_index(
            c - the_function_differential_size_);

        the_value_algebraic_[an_index] += a_velocity_sum;
    }

    the_last_time_ = a_current_time;
}

void DAESolver::initialize(Real variable_array_differential[], 
    Real variable_array_algebraic[], 
    const Integer a_variable_array_size,
    const Integer a_variable_array_size_differential)
{
    the_value_differential_.resize(a_variable_array_size_differential);
    std::copy(variable_array_differential,
              variable_array_differential + a_variable_array_size_differential,
              the_value_differential_.begin());
    the_value_differential_buffer_ = the_value_differential_;

    const VariableArray::size_type n_algebraic(
        a_variable_array_size - a_variable_array_size_differential);
    the_value_algebraic_.resize(n_algebraic);
    std::copy(variable_array_algebraic, variable_array_algebraic + n_algebraic,
              the_value_algebraic_.begin());
    the_value_algebraic_buffer_ = the_value_algebraic_;

    the_last_time_ = the_current_time_;

    is_interrupted_ = true;
    the_taylor_series_.resize(boost::extents[5]
        [static_cast< RealMatrix::index >(a_variable_array_size)]);

    eta_ = 1.0;
    the_stopping_criterion_ = std::max(10.0 * uround_ / rtoler_,
        std::min(0.03, sqrt(rtoler_)));
  
    the_first_step_flag_ = true;
    the_jacobian_calculate_flag_ = true;


    const VariableArray::size_type a_size(a_variable_array_size);
 
    if (the_system_size_ != a_size)
    {
        the_system_size_ = a_variable_array_size;
        the_function_differential_size_ = a_variable_array_size_differential;

        the_activity_algebraic_buffer_.clear();
        the_activity_algebraic_buffer_.resize(
            the_system_size_ - the_function_differential_size_);

        the_jacobian_.resize(a_size);
        for (VariableArray::size_type c(0); c < a_size; c++)
        {
            the_jacobian_[c].resize(a_size);
        }

        if (the_jacobian_matrix1_)
        {
            gsl_matrix_free(the_jacobian_matrix1_);
            the_jacobian_matrix1_ = NULLPTR;
        }
	
        if (a_size > 0)
        {
            the_jacobian_matrix1_ = gsl_matrix_calloc(a_size, a_size);
        }

        if (the_permutation1_)
        {
            gsl_permutation_free(the_permutation1_);
            the_permutation1_ = NULLPTR;
        }
	
        if (a_size > 0)
        {
            the_permutation1_ = gsl_permutation_alloc(a_size);
        }

        if (the_velocity_vector1_)
        {
            gsl_vector_free(the_velocity_vector1_);
            the_velocity_vector1_ = NULLPTR;
        }
        if (a_size > 0)
        {
            the_velocity_vector1_ = gsl_vector_calloc(a_size);
        }

        if (the_solution_vector1_)
        {
            gsl_vector_free(the_solution_vector1_);
            the_solution_vector1_ = NULLPTR;
        }
        if (a_size > 0)
        {
            the_solution_vector1_ = gsl_vector_calloc(a_size);
        }

        the_w_.resize(a_size * 3);

        if (the_jacobian_matrix2_)
        {
            gsl_matrix_complex_free(the_jacobian_matrix2_);
            the_jacobian_matrix2_ = NULLPTR;
        }
        if (a_size > 0)
        {
            the_jacobian_matrix2_ = gsl_matrix_complex_calloc(a_size, a_size);
        }

        if (the_permutation2_)
        {
            gsl_permutation_free(the_permutation2_);
            the_permutation2_ = NULLPTR;
        }
        if (a_size > 0)
        {
            the_permutation2_ = gsl_permutation_alloc(a_size);
        }

        if (the_velocity_vector2_)
        {
            gsl_vector_complex_free(the_velocity_vector2_);
            the_velocity_vector2_ = NULLPTR;
        }
        if (a_size > 0)
        {
            the_velocity_vector2_ = gsl_vector_complex_calloc(a_size);
        }

        if (the_solution_vector2_)
        {
            gsl_vector_complex_free(the_solution_vector2_);
            the_solution_vector2_ = NULLPTR;
        }
        if (a_size > 0)
        {
            the_solution_vector2_ = gsl_vector_complex_calloc(a_size);
        }
    }
}


const Real DAESolver::Interpolant::
get_difference(Time a_time, Real an_interval, 
    VariableArray::size_type the_index) const
{
    if (!the_solver_->the_state_flag_)
    {
        return 0.0;
    }

    const Real a_time_interval1(a_time - the_solver_->get_current_time());
    const Real a_time_interval2(a_time_interval1 - an_interval);

    const RealMatrix & a_taylor_series(the_solver_->get_taylor_series());
    const Real* a_taylor_coefficient_ptr(a_taylor_series.origin() + the_index);

    // calculate first order.
    // here it assumes that always a_taylor_series.size() >= 1

    // *a_taylor_coefficient_ptr := a_taylor_series[0][the_index]
    Real a_value1(*a_taylor_coefficient_ptr * a_time_interval1);
    Real a_value2(*a_taylor_coefficient_ptr * a_time_interval2);


    // check if second and higher order calculations are necessary.
    const RealMatrix::size_type a_taylor_size(3);
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
            
            // a_taylor_series[s][the_index]
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

Time DAESolver::reschedule()
{
    const Time a_local_time(the_current_time_);
    const Time a_new_step_interval(get_step_interval());

    return a_new_step_interval + a_local_time;
}

Integer DAESolver::step()
{
    return update_internal_state(the_next_step_interval_);
}


Integer DAESolver::update_internal_state(Real a_step_interval)
{
    const VariableArray::size_type a_size(the_system_size_);

    the_state_flag_ = false;
    Real const a_previous_step_interval(get_step_interval());

    clear_variables();

    the_rejected_step_flag_ = false;

    for (FunctionArray::size_type c(the_function_differential_size_);
        c < the_system_size_; c++)
    {
        const FunctionArray::size_type an_index(
            c - the_function_differential_size_);

        the_activity_algebraic_buffer_[an_index] 
            = (*the_function_[c])(the_value_differential_,
                the_value_algebraic_, the_current_time_);
    }

    set_variable_velocity(the_taylor_series_[3]);

    if (the_jacobian_calculate_flag_)
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

    UnsignedInteger count(0);
    for (;;)
    {
        std::pair< bool, Real > a_result(calculate_radauIIA(
            a_step_interval, a_previous_step_interval));

        if (a_result.first)
        {
            break;
        }

        if (++count >= 3)
        {
            break;
        }

        //comment: E-Cell3では
        //aStepInterval = a_result.first;
        //となっていたものを修正
        a_step_interval = a_result.second;

        the_rejected_step_flag_ = true;

        if (!the_jacobian_calculate_flag_)
        {
            calculate_jacobian();
            the_jacobian_calculate_flag_ = true;
        }

        set_jacobian_matrix(a_step_interval);
    }

    the_tolerable_step_interval_ = a_step_interval;


    // the_w_ will already be transformed to Z-form

    for (VariableArray::size_type c(0); c < a_size; ++c)
    {
        the_taylor_series_[3][c] = the_w_[c + a_size * 2];
        the_taylor_series_[3][c] /= a_step_interval;

        if (c < the_function_differential_size_)
        {
            the_value_differential_[c] = the_value_differential_buffer_[c];
        }
        else
        {
            const VariableArray::size_type an_index(
                c - the_function_differential_size_);
            the_value_algebraic_[an_index]
                = the_value_algebraic_buffer_[an_index];
        }
    }

    for (VariableArray::size_type c(0); c < a_size; c++)
    {
        const Real z1(the_w_[c]);
        const Real z2(the_w_[c + a_size]);
        const Real z3(the_w_[c + a_size * 2]);

        the_taylor_series_[0][c] = (13.0 + 7.0 * SQRT6) / 3.0 * z1
            + (13.0 - 7.0 * SQRT6) / 3.0 * z2
            + 1.0 / 3.0 * z3;
        the_taylor_series_[1][c] = - (23.0 + 22.0 * SQRT6) / 3.0 * z1
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

    update_internal_state_differential_stepper(a_step_interval);

    //20100224追加
    return get_status_code();

}

//20100224追加
Integer DAESolver::get_status_code()
{
    Integer status_code(0);
    Real minimum_occurrence_time(INF);
    const RealMatrix::size_type a_taylor_size(3);
    RealVector tmp_min;
    tmp_min.resize(0);

    for (StatusEventArray::size_type i(0); i < the_status_event_.size(); ++i)
    {

        Integer index = the_status_event_[i].variable_index;
        Integer flag = the_status_event_[i].variable_flag;
        Real threshold = the_status_event_[i].threshold;
        Integer code = the_status_event_[i].status_code;

        Real a[a_taylor_size + 1];

        if (flag == 0)
        {
            a[0] = the_value_differential_[index] - threshold;
            a[1] = the_taylor_series_[0][index];
            Real power = 1.0;
            for (Integer j(0); j < a_taylor_size - 1; ++j)
            {
                power = power * the_tolerable_step_interval_;
                a[j + 2] = the_taylor_series_[j + 1][index] / power;
            }
        }
        else
        {
            a[0] = the_value_algebraic_[index] - threshold;
            a[1] = the_taylor_series_[0][
                index + the_function_differential_size_];
            Real power = 1.0;
            for (Integer j(0); j < a_taylor_size - 1; ++j)
            {
                power = power * the_tolerable_step_interval_;
                a[j + 2] = the_taylor_series_[j + 1][
                    index + the_function_differential_size_] / power;
            }
        }

        Real z[2 * a_taylor_size];
        gsl_poly_complex_workspace * w
            = gsl_poly_complex_workspace_alloc(a_taylor_size + 1);
        gsl_poly_complex_solve(a, a_taylor_size + 1, w, z);
        gsl_poly_complex_workspace_free(w);

        for (Integer k(0); k < a_taylor_size; k++)
        {
            const Real solution(z[2 * k]);

            if (z[2 * k + 1] == 0.0 && solution >= 0. &&
                solution < the_next_time_ - the_current_time_)
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
        if (tmp_min[n] == minimum_occurrence_time) count++;
    }
    if (count > 1)
    {
        throw "Some status events occurred at the same time.";
    }

    return status_code;
}

void DAESolver::clear_variables()
{
    // save original value values
    the_value_differential_buffer_.clear();
    the_value_differential_buffer_ = the_value_differential_;
    the_value_algebraic_buffer_.clear();
    the_value_algebraic_buffer_ = the_value_algebraic_;
}



void DAESolver::set_variable_velocity(
    boost::detail::multi_array::sub_array<Real, 1> a_velocity_buffer)
{
    for (RealMatrix::index i(0);
        i < static_cast< RealMatrix::index >(a_velocity_buffer.size());
        ++i)
    {
        a_velocity_buffer[i] = 0.0;
    }

    for (FunctionArray::size_type i(0);
        i < the_function_differential_size_; ++i)
    {
        a_velocity_buffer[i] = (*the_function_[i])(
            the_value_differential_, the_value_algebraic_, the_current_time_);
    }

}


void DAESolver::calculate_jacobian()
{
    Real a_perturbation;

    for (VariableArray::size_type i(0); i < the_system_size_; ++i)
    {
        Real a_tmp_value(0.0);
        if (i < the_function_differential_size_)
        {
            a_tmp_value = the_value_differential_[i];
        }
        else
        {
            const VariableArray::size_type an_index(
                i - the_function_differential_size_);
            a_tmp_value = the_value_algebraic_[an_index];
        }
        const Real a_value(a_tmp_value);

        a_perturbation = sqrt(uround_ * std::max(1e-5, fabs(a_value)));

        if (i < the_function_differential_size_)
        {
            the_value_differential_[i]
                = the_value_differential_buffer_[i] + a_perturbation;
        }
        else
        {
            const VariableArray::size_type an_index(
                i - the_function_differential_size_);
            the_value_algebraic_[an_index]
                = the_value_algebraic_buffer_[an_index] + a_perturbation;
        } 

        set_variable_velocity(the_taylor_series_[4]);

        for (FunctionArray::size_type c(the_function_differential_size_);
            c < the_function_.size(); c++)
        {
            const FunctionArray::size_type an_index(
                c - the_function_differential_size_);

            the_jacobian_[c][i] = -((*the_function_[c])(
                the_value_differential_, the_value_algebraic_,
                the_current_time_)
                - the_activity_algebraic_buffer_[an_index]) /a_perturbation;
        }

        for (VariableArray::size_type j(0);
            j < the_function_differential_size_; ++j)
        {
            // int as VariableVector::size_type
            the_jacobian_[j][i] = - (the_taylor_series_[4][j]
                - the_taylor_series_[3][j]) / a_perturbation;
        }

        if (i < the_function_differential_size_)
        {
            the_value_differential_[i] = a_value;
        }
        else
        {
            const VariableArray::size_type an_index(
                i - the_function_differential_size_);
            the_value_algebraic_[an_index] = a_value;
        }

    }

}


void DAESolver::set_jacobian_matrix(Real const a_step_interval)
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

    for (Integer c(0); c < the_function_differential_size_; ++c)
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

void DAESolver::decomp_jacobian_matrix()
{
    int a_sign_num;

    if (the_system_size_ == 0)
    {
        return;
    }

    gsl_linalg_LU_decomp(the_jacobian_matrix1_,
        the_permutation1_, &a_sign_num);
    gsl_linalg_complex_LU_decomp(the_jacobian_matrix2_,
        the_permutation2_, &a_sign_num);
}

std::pair< bool, Real > DAESolver::calculate_radauIIA(
    Real const a_step_interval, Real const a_previous_step_interval)
{
    const VariableArray::size_type a_size(the_system_size_);
    Real a_new_step_interval;
    Real a_norm;
    Real theta(fabs(the_jacobian_recalculate_theta_));

    Integer an_iterator(0);

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
            const Real cont1(the_taylor_series_[0][c]
                + 2.0 * cont2 - 3.0 * cont3);

            const Real z1(a_previous_step_interval
                * c1q * (cont1 + c1q * (cont2 + c1q * cont3)));
            const Real z2(a_previous_step_interval
                * c2q * (cont1 + c2q * (cont2 + c2q * cont3)));
            const Real z3(a_previous_step_interval
                * c3q * (cont1 + c3q * (cont2 + c3q * cont3)));

            the_w_[c] = 4.3255798900631553510 * z1
                + 0.33919925181580986954 * z2 + 0.54177053993587487119 * z3;
            the_w_[c + a_size] = -4.1787185915519047273 * z1
                - 0.32768282076106238708 * z2 + 0.47662355450055045196 * z3;
            the_w_[c + a_size * 2] = -0.50287263494578687595 * z1
                + 2.5719269498556054292 * z2 - 0.59603920482822492497 * z3;
        }
    }
    else
    {
        for (VariableArray::size_type c(0); c < the_system_size_; ++c)
        {
            the_w_[c] = 0.0;
            the_w_[c + the_system_size_] = 0.0;
            the_w_[c + the_system_size_ * 2] = 0.0;
        }
    }

    eta_ = pow(std::max(eta_, uround_), 0.8);

    while (an_iterator < the_max_iteration_number_)
    {
        calculate_rhs(a_step_interval);

        const Real a_previous_norm(std::max(a_norm, uround_));
        a_norm = solve();

        if (an_iterator > 0 && (an_iterator != the_max_iteration_number_ - 1))
        {
            const Real a_theta_q = a_norm / a_previous_norm;
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
                const Real an_iteration_error(eta_ * 
                    a_norm * 
                    pow(theta, static_cast<int>(
                    the_max_iteration_number_ - 2 - an_iterator)) /
                    the_stopping_criterion_);

                if (an_iteration_error >= 1.0)
                {
                    return std::make_pair(false, a_step_interval * 0.8 *
                        std::pow(std::max(1e-4, std::min(20.0,
                        an_iteration_error)),
                        - 1.0 / (4 + the_max_iteration_number_
                        - 2 - an_iterator)));
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
    for (VariableArray::size_type c(0); c < a_size; ++c)
    {
        const Real w1(the_w_[c]);
        const Real w2(the_w_[c + a_size]);
        const Real w3(the_w_[c + a_size * 2]);

        the_w_[c] = w1 * 0.091232394870892942792
            - w2 * 0.14125529502095420843
            - w3 * 0.030029194105147424492;
        the_w_[c + a_size] = w1 * 0.24171793270710701896
            + w2 * 0.20412935229379993199
            + w3 * 0.38294211275726193779;
        the_w_[c + a_size * 2] = w1 * 0.96604818261509293619 + w2;
    }

    const Real an_error(estimate_local_error(a_step_interval));

    Real a_safety_factor(std::min(0.9, 0.9 *
        (1 + 2 * the_max_iteration_number_) /
        (an_iterator + 1 + 2 * the_max_iteration_number_)));
    a_safety_factor = std::max(0.125,
        std::min(5.0, pow(an_error, 0.25) / a_safety_factor));

    a_new_step_interval = a_step_interval / a_safety_factor;

    if (an_error < 1.0)
    {
        // step is accepted
        if (!the_first_step_flag_)
        {
            Real a_gustafsson_factor(
                the_accepted_step_interval_ / a_step_interval
                * pow(an_error * an_error / the_accepted_error_, 0.25) / 0.9);
            a_gustafsson_factor
                = std::max(0.125, std::min(5.0, a_gustafsson_factor));

            if (a_safety_factor < a_gustafsson_factor)
            {
                a_safety_factor = a_gustafsson_factor;
                a_new_step_interval = a_step_interval / a_gustafsson_factor;
            }
        }

        the_accepted_step_interval_ = a_step_interval;
        the_accepted_error_ = std::max(1.0e-2, an_error);

        if (the_rejected_step_flag_)
        {
            a_new_step_interval
                = std::min(a_new_step_interval, a_step_interval);
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
            the_next_step_interval_
                = std::min(a_step_interval, the_max_step_interval_);
        }
        else
        {
            the_next_step_interval_
                = std::min(a_new_step_interval, the_max_step_interval_);
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


void DAESolver::calculate_rhs(Real a_step_interval)
{
    const Time a_current_time(the_current_time_);
    const VariableArray::size_type a_size(the_system_size_);

    const Real alphah(alpha_ / a_step_interval);
    const Real betah(beta_ / a_step_interval);
    const Real gammah(gamma_ / a_step_interval);

    gsl_complex comp;

    RealVector aTIF;
    aTIF.resize(the_system_size_ * 3);

    for (VariableArray::size_type c(0); c < the_system_size_; ++c)
    {
        const Real z(the_w_[c] * 0.091232394870892942792
            - the_w_[c + a_size] * 0.14125529502095420843
            - the_w_[c + 2 * a_size] * 0.030029194105147424492);

        if (c < the_function_differential_size_)
        the_value_differential_[c] = the_value_differential_buffer_[c] + z;
        else
        {
            const VariableArray::size_type an_index(
                c - the_function_differential_size_);
            the_value_algebraic_[an_index]
                = the_value_algebraic_buffer_[an_index] + z;
        }
    }

    // ========= 1 ===========
    the_current_time_
        = a_current_time + a_step_interval * (4.0 - SQRT6) / 10.0;
    set_variable_velocity(the_taylor_series_[4]);

    for (FunctionArray::size_type c(
        the_function_differential_size_); c < the_system_size_; c++)
    {
        Real tmp_value = (*the_function_[c])(the_value_differential_,
            the_value_algebraic_, the_current_time_);
        aTIF[c] = tmp_value * 4.3255798900631553510;
        aTIF[c + a_size] = tmp_value * -4.1787185915519047273;
        aTIF[c + a_size * 2] = tmp_value * -0.50287263494578687595;
    }

    for (VariableArray::size_type c(0);
        c < the_function_differential_size_; c++)
    {
        aTIF[c] = the_taylor_series_[4][c] * 4.3255798900631553510;
        aTIF[c + a_size] = the_taylor_series_[4][c] * -4.1787185915519047273;
        aTIF[c + a_size * 2]
            = the_taylor_series_[4][c] * -0.50287263494578687595;
    }

    for (VariableArray::size_type c(0); c < a_size; ++c)
    {
        const Real z(the_w_[c] * 0.24171793270710701896
            + the_w_[c + a_size] * 0.20412935229379993199
            + the_w_[c + 2 * a_size] * 0.38294211275726193779);

        if (c < the_function_differential_size_)
        {
            the_value_differential_[c]
                = the_value_differential_buffer_[c] + z;
        }
        else
        {
            const VariableArray::size_type an_index(
                c - the_function_differential_size_);
            the_value_algebraic_[an_index]
                = the_value_algebraic_buffer_[an_index] + z;
        }
    }

    // ========= 2 ===========

    the_current_time_
        = a_current_time + a_step_interval * (4.0 + SQRT6) / 10.0;
    set_variable_velocity(the_taylor_series_[4]);

    for (FunctionArray::size_type c(the_function_differential_size_);
        c < the_system_size_; c++)
    {
        Real tmp_value = (*the_function_[c])(
            the_value_differential_, the_value_algebraic_, the_current_time_);
        aTIF[c] += tmp_value * 0.33919925181580986954;
        aTIF[c + a_size] -= tmp_value * 0.32768282076106238708;
        aTIF[c + a_size * 2] += tmp_value * 2.5719269498556054292;
    }

    for (VariableArray::size_type c(0);
        c < the_function_differential_size_; c++)
    {
        aTIF[c] += the_taylor_series_[4][c] * 0.33919925181580986954;
        aTIF[c + a_size] -= the_taylor_series_[4][c] * 0.32768282076106238708;
        aTIF[c + a_size * 2]
            += the_taylor_series_[4][c] * 2.5719269498556054292;
    }

    for (VariableArray::size_type c(0); c < a_size; ++c)
    {
        const Real z(the_w_[c] * 0.96604818261509293619 + the_w_[c + a_size]);

        if (c < the_function_differential_size_)
        {
            the_value_differential_[c]
                = the_value_differential_buffer_[c] + z;
        }
        else
        {
            const VariableArray::size_type an_index(
                c - the_function_differential_size_);
            the_value_algebraic_[an_index]
                = the_value_algebraic_buffer_[an_index] + z; 
        }
    }

    // ========= 3 ===========

    the_current_time_ = a_current_time + a_step_interval;
    set_variable_velocity(the_taylor_series_[4]);

    for (FunctionArray::size_type c(the_function_differential_size_);
        c < the_system_size_; c++)
    {
        Real tmp_value = (*the_function_[c])(
            the_value_differential_, the_value_algebraic_, the_current_time_);
        aTIF[c] += tmp_value * 0.54177053993587487119;
        aTIF[c + a_size] += tmp_value * 0.47662355450055045196;
        aTIF[c + a_size * 2] -= tmp_value * 0.59603920482822492497;

        gsl_vector_set(the_velocity_vector1_, c, aTIF[c]);

        GSL_SET_COMPLEX(&comp, aTIF[c + a_size], aTIF[c + a_size * 2]);
        gsl_vector_complex_set(the_velocity_vector2_, c, comp);
    }

    for (VariableArray::size_type c(0);
        c < the_function_differential_size_; c++)
    {
        aTIF[c] += the_taylor_series_[4][c] * 0.54177053993587487119;
        aTIF[c + a_size] += the_taylor_series_[4][c] * 0.47662355450055045196;
        aTIF[c + a_size * 2]
            -= the_taylor_series_[4][c] * 0.59603920482822492497;

        gsl_vector_set(the_velocity_vector1_,
            c, aTIF[c] - the_w_[c] * gammah);

        GSL_SET_COMPLEX(&comp,
            aTIF[c + a_size]
            - the_w_[c + a_size] * alphah
            + the_w_[c + a_size * 2] * betah,
            aTIF[c + a_size * 2]
            - the_w_[c + a_size] * betah
            - the_w_[c + a_size * 2] * alphah);
        gsl_vector_complex_set(the_velocity_vector2_, c, comp);
    }

    the_current_time_ = a_current_time;
}


Real DAESolver::solve()
{
    if (the_system_size_ == 0)
    {
        return 0.0;
    }

    const VariableArray::size_type a_size(the_system_size_);

    gsl_linalg_LU_solve(the_jacobian_matrix1_, the_permutation1_,
        the_velocity_vector1_, the_solution_vector1_);
    gsl_linalg_complex_LU_solve(the_jacobian_matrix2_, the_permutation2_,
        the_velocity_vector2_, the_solution_vector2_);

    Real a_norm(0.0);
    Real delta_w(0.0);
    gsl_complex comp;

    VariableArray a_tmp_value_buffer = the_value_differential_buffer_;
    a_tmp_value_buffer.insert(
        a_tmp_value_buffer.end(), the_value_algebraic_buffer_.begin(),
        the_value_algebraic_buffer_.end());

    for (VariableArray::size_type c(0); c < a_size; ++c)
    {
        Real a_tolerance2(rtoler_ * fabs(a_tmp_value_buffer[c]) + atoler_);

        a_tolerance2 = a_tolerance2 * a_tolerance2;

        delta_w = gsl_vector_get(the_solution_vector1_, c);
        the_w_[c] += delta_w;
        a_norm += delta_w * delta_w / a_tolerance2;

        comp = gsl_vector_complex_get(the_solution_vector2_, c);

        delta_w = GSL_REAL(comp);
        the_w_[c + a_size] += delta_w;
        a_norm += delta_w * delta_w / a_tolerance2;

        delta_w = GSL_IMAG(comp);
        the_w_[c + a_size * 2] += delta_w;
        a_norm += delta_w * delta_w / a_tolerance2;
    }

    return sqrt(a_norm / (3 * a_size));
}

Real DAESolver::estimate_local_error(Real const a_step_interval)
{
    if (the_system_size_ == 0)
    {
        return 0.0;
    }

    const VariableArray::size_type a_size(the_system_size_);

    Real an_error;

    const Real hee1((-13.0 - 7.0 * SQRT6) / (3.0 * a_step_interval));
    const Real hee2((-13.0 + 7.0 * SQRT6) / (3.0 * a_step_interval));
    const Real hee3(-1.0 / (3.0 * a_step_interval));

    // the_w_ will already be transformed to Z-form
    for (VariableArray::size_type c(the_function_differential_size_);
        c < the_system_size_; c++)
    {
        const VariableArray::size_type an_index(
            c - the_function_differential_size_);
        gsl_vector_set(the_velocity_vector1_, c,
            the_activity_algebraic_buffer_[an_index]);
    }

    for (VariableArray::size_type c(0);
        c < the_function_differential_size_; c++)
    {
        gsl_vector_set(the_velocity_vector1_, c,
            the_taylor_series_[3][c]
            + the_w_[c] * hee1
            + the_w_[c + a_size] * hee2
            + the_w_[c + 2 * a_size] * hee3);
    }

    gsl_linalg_LU_solve(the_jacobian_matrix1_, the_permutation1_,
            the_velocity_vector1_, the_solution_vector1_);

    an_error = 0.0;

    VariableArray a_tmp_value_buffer = the_value_differential_buffer_;
    a_tmp_value_buffer.insert(a_tmp_value_buffer.end(),
        the_value_algebraic_buffer_.begin(),
        the_value_algebraic_buffer_.end());

    for (VariableArray::size_type c(0); c < a_size; ++c)
    {
        const Real a_tolerance(rtoler_
            * fabs(a_tmp_value_buffer[c]) + atoler_);

        Real a_difference(gsl_vector_get(the_solution_vector1_, c));

        if (c < the_function_differential_size_)
        {
            the_value_differential_[c]
                = the_value_differential_buffer_[c] + a_difference;
        }
        else
        {
            const VariableArray::size_type an_index(
                c - the_function_differential_size_);
            the_value_algebraic_[an_index]
                = the_value_algebraic_buffer_[an_index] + a_difference;
        }

        a_difference /= a_tolerance;
        an_error += a_difference * a_difference;
    }

    an_error = std::max(sqrt(an_error / a_size), 1e-10);

    if (an_error < 1.0)
    {
        return an_error;
    }

    if (the_first_step_flag_ || the_rejected_step_flag_)
    {
        set_variable_velocity(the_taylor_series_[4]);

        for (FunctionArray::size_type c(the_function_differential_size_);
            c < the_system_size_; c++)
        {
            gsl_vector_set(the_velocity_vector1_, c,
                (*the_function_[c])(the_value_differential_,
                the_value_algebraic_, the_current_time_));
        }

        for (VariableArray::size_type c(0);
            c < the_function_differential_size_; c++)
        {
            gsl_vector_set(the_velocity_vector1_, c,
                the_taylor_series_[4][c]
                + the_w_[c] * hee1
                + the_w_[c + a_size] * hee2
                + the_w_[c + 2 * a_size] * hee3);
        }

        gsl_linalg_LU_solve(the_jacobian_matrix1_, the_permutation1_,
            the_velocity_vector1_, the_solution_vector1_);

        an_error = 0.0;
        for (VariableArray::size_type c(0); c < a_size; ++c)
        {
            const Real a_tolerance(rtoler_ * fabs(a_tmp_value_buffer[c])
                + atoler_);

            Real a_difference(gsl_vector_get(the_solution_vector1_, c));
            a_difference /= a_tolerance;

            an_error += a_difference * a_difference;
        }

        an_error = std::max(sqrt(an_error / a_size), 1e-10);

    }

    return an_error;
}

void DAESolver::update_internal_state_differential_stepper(
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
