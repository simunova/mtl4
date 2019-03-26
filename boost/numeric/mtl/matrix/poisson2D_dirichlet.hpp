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

#ifndef MTL_MATRIX_POISSON2D_DIRICHLET_INCLUDE
#define MTL_MATRIX_POISSON2D_DIRICHLET_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/vector/mat_cvec_multiplier.hpp>

namespace mtl { namespace mat {

/// Matrix-free linear operator for a Poisson equation on a rectangular domain of \p m by \p n with Dirichlet boundary conditions
struct poisson2D_dirichlet
{
    /// Constructor
    poisson2D_dirichlet(int m, int n) : m(m), n(n), s(m * n) {}

    /// Member function that realizes the multiplication
    template <typename VectorIn, typename VectorOut, typename Assign>
    void mult(const VectorIn& v, VectorOut& w, Assign) const
    {
	MTL_DEBUG_THROW_IF(int(size(v)) != s, incompatible_size());
	MTL_DEBUG_THROW_IF(size(v) != size(w), incompatible_size());

	const int nb = n < 3 ? 1 : (n - 2) / 4 * 4 + 1;

	// Inner domain
	for (int i= 1; i < m-1; i++) {
	    int kmax= i * n + nb;
	    for (int k= i * n + 1; k < kmax; k+= 4) {
		typename Collection<VectorIn>::value_type const v0= v[k], v1= v[k+1], v2= v[k+2], v3= v[k+3];
		Assign::apply(w[k], 4 * v0 - v[k-n] - v[k+n] - v[k-1] - v1); 
		Assign::apply(w[k+1], 4 * v1 - v[k-n+1] - v[k+n+1] - v0 - v2); 
		Assign::apply(w[k+2], 4 * v2 - v[k-n+2] - v[k+n+2] - v1 - v3); 
		Assign::apply(w[k+3], 4 * v3 - v[k-n+3] - v[k+n+3] - v2 - v[k+4]); 
	    }
	    for (int j= nb, k= i * n + j; j < n-1; j++, k++) 
		Assign::apply(w[k], 4 * v[k] - v[k-n] - v[k+n] - v[k-1] - v[k+1]); 
	}
	    
	// Upper border
	for (int j= 1; j < n-1; j++) 
	    Assign::apply(w[j], 4 * v[j] - v[j+n] - v[j-1] - v[j+1]);

	// Lower border
	for (int j= 1, k= (m-1) * n + j; j < n-1; j++, k++) 
	    Assign::apply(w[k], 4 * v[k] - v[k-n] - v[k-1] - v[k+1]); 
	
	// Left border
	for (int i= 1, k= n; i < m-1; i++, k+= n)
	    Assign::apply(w[k], 4 * v[k] - v[k-n] - v[k+n] - v[k+1]); 

	// Right border
	for (int i= 1, k= n+n-1; i < m-1; i++, k+= n)
	    Assign::apply(w[k], 4 * v[k] - v[k-n] - v[k+n] - v[k-1]); 

	// Corners
	Assign::apply(w[0], 4 * v[0] - v[1] - v[n]);
	Assign::apply(w[n-1], 4 * v[n-1] - v[n-2] - v[2*n - 1]);
	Assign::apply(w[(m-1)*n], 4 * v[(m-1)*n] - v[(m-2)*n] - v[(m-1)*n+1]);
	Assign::apply(w[m*n-1], 4 * v[m*n-1] - v[m*n-2] - v[m*n-n-1]);
    }

    /// Multiplication is procastinated until we know where the product goes
    template <typename VectorIn>
    vec::mat_cvec_multiplier<poisson2D_dirichlet, VectorIn> operator*(const VectorIn& v) const
    {	return vec::mat_cvec_multiplier<poisson2D_dirichlet, VectorIn>(*this, v);    }

    int m, n, s;
};

inline std::size_t size(const poisson2D_dirichlet& A) { return A.s * A.s; } ///< Matrix size
inline std::size_t num_rows(const poisson2D_dirichlet& A) { return A.s; } ///< Number of rows
inline std::size_t num_cols(const poisson2D_dirichlet& A) { return A.s; } ///< Number of columns

}} // namespace mtl::matrix

namespace mtl { 

    template <>
    struct Collection<mat::poisson2D_dirichlet>
    {
	typedef double value_type;
	typedef int    size_type;
    };

    namespace ashape {
	template <> struct ashape_aux<mtl::mat::poisson2D_dirichlet> 
	{	typedef nonscal type;    };
    }
}

#endif // MTL_MATRIX_POISSON2D_DIRICHLET_INCLUDE
