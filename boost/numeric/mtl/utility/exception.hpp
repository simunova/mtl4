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

#ifndef MTL_MTL_EXCEPTION_INCLUDE
#define MTL_MTL_EXCEPTION_INCLUDE

#include <cassert>
#include <stdexcept>

namespace mtl {

#ifndef NDEBUG
#  define MTL_DEBUG_ARG(Arg) Arg
#else
#  define MTL_DEBUG_ARG(Arg)
#endif

#ifndef MTL_ASSERT_FOR_THROW
#  define MTL_THROW_ARG(Arg) Arg
#else
#  define MTL_THROW_ARG(Arg)
#endif

// If MTL_ASSERT_FOR_THROW is defined all throws become assert
// MTL_DEBUG_THROW_IF completely disappears if NDEBUG is defined
#ifndef NDEBUG
#  ifdef MTL_ASSERT_FOR_THROW
#    define MTL_DEBUG_THROW_IF(Test, Exception) \
        { assert(!(Test)); }
#  else
#    define MTL_DEBUG_THROW_IF(Test, Exception) \
        { if (Test) throw Exception; }
#  endif
#else
#  define MTL_DEBUG_THROW_IF(Test,Exception)
#endif


#if defined(MTL_ASSERT_FOR_THROW) && !defined(NDEBUG)
#  define MTL_THROW_IF(Test, Exception)       \
   {                                          \
       assert(!(Test));			      \
   }
#else
#  define MTL_THROW_IF(Test, Exception)       \
   {                                          \
      if (Test) throw Exception;              \
   }
#endif


#if defined(MTL_ASSERT_FOR_THROW) && !defined(NDEBUG)
#  define MTL_THROW(Exception)       \
   {                                 \
       assert(0);		     \
   }
#else
#  define MTL_THROW(Exception)       \
   {                                 \
      throw Exception;               \
   }
#endif


#if 0 
standard errors:

exception
    logic_error
        domain_error
        invalid_argument
        length_error
        out_of_range
    runtime_error
        range_error
        overflow_error
        underflow_error
bad_alloc
bad_cast
bad_exception
bad_typeid

#endif

/// Exception for indices out of range
struct index_out_of_range : public std::out_of_range
{
    /// Error can be specified more precisely in constructor if desired
    explicit index_out_of_range(const char *s= "Index out of range") : std::out_of_range(s) {}
};

/// Exception for invalid range definitions, esp. in constructors
struct range_error : public std::range_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit range_error(const char *s= "Invalid range") : std::range_error(s) {}
};

/// Domain errors in MTL4
struct domain_error : public std::domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit domain_error(const char *s= "MTL4 domain error.") : std::domain_error(s) {}
};

/// Exception for arguments with incompatible sizes
struct incompatible_size : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit incompatible_size(const char *s= "Arguments have incompatible size.")
	: domain_error(s) {}
};

/// Exception for arguments that shall not be empty 
struct need_nonempty : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit need_nonempty(const char *s= "Argument must be non-empty.")
	: domain_error(s) {}
};

/// Exception for trying to change a fixed size (to another value)
struct change_static_size : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit change_static_size(const char *s= "You try to change a fixed size (to another value).")
	: domain_error(s) {}
};

/// Exception for arguments with incompatible shapes, e.g. adding matrices and vectors
struct argument_result_conflict : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit argument_result_conflict(const char *s= "Used same object illegally as argument and result.")
	: domain_error(s) {}
};

/// Exception for arguments with incompatible shapes, e.g. adding matrices and vectors
struct incompatible_shape : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit incompatible_shape(const char *s= "Arguments have incompatible shape.")
	: domain_error(s) {}
};

/// Exception for arguments with incompatible sizes
struct matrix_not_square : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit matrix_not_square(const char *s= "Matrix must be square for this operation.")
	: domain_error(s) {}
};

/// Exception for matrices too small for certain algorithms
struct matrix_too_small : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit matrix_too_small(const char *s= "Matrix is too small for certain algorithms.")
	: domain_error(s) {}
};

/// Exception for singular matrices in solvers
struct matrix_singular : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit matrix_singular(const char *s= "Matrix singular in solver.")
	: domain_error(s) {}
};

/// Exception for arguments with incompatible sizes
struct missing_diagonal : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit missing_diagonal(const char *s= "Diagonal entry missing or not where it belongs to.")
	: domain_error(s) {}
};

/// Accessing (illegally) matrix or vector during insertion phase (dense non-distributed can be accessed always)
struct access_during_insertion : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit access_during_insertion(const char *s= "Diagonal entry missing.")
	: domain_error(s) {}
};

/// Exception for a result that is not what it should be
struct unexpected_result : public domain_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit unexpected_result(const char *s= "The result of an operation is not the expected one.")
	: domain_error(s) {}
};

/// Exception for run-time errors that doesn't fit into specific categories
struct runtime_error : public std::runtime_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit runtime_error(const char *s= "Run-time error") : std::runtime_error(s) {}
};

struct division_by_zero : mtl::runtime_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit division_by_zero(const char *s= "Division by zero") : runtime_error(s) {}
};


/// Exception for logic errors that doesn't fit into specific categories
struct logic_error : public std::logic_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit logic_error(const char *s= "Logic error") : std::logic_error(s) {}
};

/// Exception for I/O errors in general
struct io_error : public std::runtime_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit io_error(const char *s= "I/O error") : std::runtime_error(s) {}
};

/// File not found
struct file_not_found : public io_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit file_not_found(const char *s= "File not found") : io_error(s) {}
};




} // namespace mtl

#endif // MTL_MTL_EXCEPTION_INCLUDE
