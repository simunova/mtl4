// Software License for MTL
//
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
//
// This file is part of the Matrix Template Library
//
// See also license.mtl.txt in the distribution.

#ifndef ITL_BASIC_ITERATION_INCLUDE
#define ITL_BASIC_ITERATION_INCLUDE

#include <iostream>
#include <complex>
#include <string>

namespace itl {

/// Basic utility class to control iterative solvers
template <class Real>
class basic_iteration
{
  public:
    typedef basic_iteration self;
    typedef Real            real;

    /// Constructor
    template <class Vector>
    basic_iteration(const Vector& r0, int max_iter_, Real t, Real a = Real(0))
      : error(0), i(0), my_norm_r0(abs(two_norm(r0))),
	max_iter(max_iter_), rtol_(t), atol_(a), is_finished(false), my_quite(false), my_suppress(false) { }

    /// Constructor
    basic_iteration(Real nb, int max_iter_, Real t, Real a = Real(0))
      : error(0), i(0), my_norm_r0(nb), max_iter(max_iter_), rtol_(t), atol_(a), is_finished(false), 
	my_quite(false), my_suppress(false) {}

    virtual ~basic_iteration() {}

    bool check_max()
    {
	if (i >= max_iter) 
	    error= 1, is_finished= true, err_msg= "Too many iterations.";
	return is_finished;
    }

    /// Iteration finished according to the norm of r
    template <class Vector>
    bool finished(const Vector& r) 
    {
	if (converged(two_norm(r)))
	    return is_finished= true;
	return check_max();
    }

    /// Iteration finished according to residual value r
    bool finished(const Real& r) 
    {
	if (converged(r))
	    return is_finished= true;
	return check_max();
    }

    /// Iteration finished according to complex residual value r
    template <typename T>
    bool finished(const std::complex<T>& r) 
    {
	if (converged(std::abs(r))) 
	    return is_finished= true;
	return check_max();
    }

    /// Iteration finished according to last provided residual
    bool finished() const { return is_finished; }

    template <class T>
    int terminate(const T& r) { finished(r); return error; }

  private:
    bool converged(const Real& r) { resid_= r; return converged(); } 

    bool converged() const 
    { 
	// std::cout << "abs: " << resid_ << " <= " << atol_ << '\n';
	// std::cout << "rel: " << resid_ << " <= " << rtol_ << " * " << my_norm_r0 
	// 	  << " = " << rtol_ * my_norm_r0 << '\n';
	if (my_norm_r0 == 0) 
	    return resid_ <= atol_;  // ignore relative tolerance if |r0| is zero
	return resid_ <= rtol_ * my_norm_r0 || resid_ <= atol_;
    }

  public:
    self& operator++() { ++i; return *this; } ///< Increment counter

    self& operator+=(int n) { i+= n; return *this; } ///< Increment counter by n

    bool first() const { return i <= 1; } ///< Is first iteration?

    virtual operator int() const { return error; } ///< Conversion to int returns error code

    virtual int error_code() const { return error; } ///< Error code

    bool is_converged() const { return is_finished && error == 0; } ///< Reached tolerance without error

    int iterations() const { return i; } ///< Already performed iterations
    
    int max_iterations() const { return max_iter; } ///< Maximal number of iterations

    void set_max_iterations(int m) { max_iter= m; } ///< Set maximal number of iterations

    void restart() { i= 0; } ///< Reset the iteration number

    Real resid() const { return resid_; } ///< Last residuum

    Real relresid() const { return resid_ / my_norm_r0; } ///< Last residuum compared with initial one

    Real normb() const { return my_norm_r0; } // deprecated

    Real norm_r0() const { return my_norm_r0; } ///< Initial residual
    void set_norm_r0(Real r0) { my_norm_r0= r0; } ///< Reset initial residual

    Real tol() const { return rtol_; } ///< Relative tolerance
    Real atol() const { return atol_; } ///< Absolute tolerance

    /// Fail with error code
    int fail(int err_code) { error = err_code; return error_code(); }

    /// Fail with error code and message
    int fail(int err_code, const std::string& msg)
    { error = err_code; err_msg = msg; return error_code(); }

    void set(Real v) { my_norm_r0 = v; } ///< Reset initial residual

    void set_quite(bool q) { my_quite= q; } ///< Turn logging off (or on)

    bool is_quite() const { return my_quite; } ///< Is logging turned off

    void suppress_resume(bool s) { my_suppress= s; } ///< Suppress final resume

    bool resume_suppressed() const { return my_suppress; }///< Is final resume suppressed

    void update_progress(const basic_iteration& that)
    {
	i= that.i;
	resid_= that.resid_;
	if (that.error > 1) { // copy error except too many iterations
	    error= that.error;
	    err_msg= that.err_msg;
	    is_finished= true;
	} else 
	    finished(resid_);
    }

  protected:
    int          error, i;
    Real         my_norm_r0;
    int          max_iter;
    Real         rtol_, atol_, resid_;
    std::string  err_msg;
    bool         is_finished, my_quite, my_suppress;
};


} // namespace itl

#endif // ITL_BASIC_ITERATION_INCLUDE
