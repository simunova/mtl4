// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_MATRIX_IDENTITY2D_INCLUDE
#define MTL_MATRIX_IDENTITY2D_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/vector/mat_cvec_multiplier.hpp>

namespace mtl { namespace mat {

/// Matrix-free linear operator for identity
struct identity2D
{
    /// Constructor for \p m by \p m matrix
    identity2D(std::size_t m) : m(m), n(m) {}

    /// Constructor for \p m by \p n matrix 
    identity2D(std::size_t m, std::size_t n) : m(m), n(n) {}

    /// Member function that realizes the multiplication
    template <typename VectorIn, typename VectorOut, typename Assign>
    void mult(const VectorIn& v, VectorOut& w, Assign) const
    {
	MTL_DEBUG_THROW_IF(std::size_t(size(v)) != n, incompatible_size());
	MTL_DEBUG_THROW_IF(size(w) != 0 && std::size_t(size(w)) != m, incompatible_size());

	if (size(w) == 0)
	    w.change_dim(m);

	if (m == n) 
	    Assign::first_update(w, v);
	else if (m < n)
	    Assign::first_update(w, v[irange(m)]);
	else {
	    VectorOut w1(w[irange(n)]), w2(w[irange(n, imax)]);
	    Assign::first_update(w1, v);
	    Assign::init(w2);
	}
    }

    /// Multiplication is procastinated until we know where the product goes
    template <typename VectorIn>
    vec::mat_cvec_multiplier<identity2D, VectorIn> operator*(const VectorIn& v) const
    {	return vec::mat_cvec_multiplier<identity2D, VectorIn>(*this, v);    }

    std::size_t m, n;
};

inline std::size_t size(const identity2D& A) { return A.m * A.n; } ///< Matrix size
inline std::size_t num_rows(const identity2D& A) { return A.m; } ///< Number of rows
inline std::size_t num_cols(const identity2D& A) { return A.n; } ///< Number of columns

}} // namespace mtl::matrix

namespace mtl { 

    template <>
    struct Collection<mat::identity2D>
    {
	typedef double         value_type;
	typedef std::size_t    size_type;
    };

    // namespace ashape {
    // 	template <> struct ashape_aux<mtl::mat::identity2D> 
    // 	{	typedef nonscal type;    };
    // }
}

#endif // MTL_MATRIX_IDENTITY2D_INCLUDE
