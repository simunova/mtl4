/*
 * Marcel Schiffel, 13.10.11
 * 
 * definition of default relaxation parameter for iterative solvers
 */


#ifndef MTL_ITL_RELAXATION_PARAMETER_HPP
#define MTL_ITL_RELAXATION_PARAMETER_HPP


namespace itl {


/**
 * default relaxation parameter for iterative solvers
 */
struct default_omega
{
    static const double value;
};

const double default_omega::value= 2./3.;

} /* namespace itl */

#endif

