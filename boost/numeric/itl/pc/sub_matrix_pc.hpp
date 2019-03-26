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

#ifndef ITL_PC_SUB_MATRIX_PC_INCLUDE
#define ITL_PC_SUB_MATRIX_PC_INCLUDE

#include <boost/static_assert.hpp>
#include <boost/mpl/if.hpp>

#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/itl/pc/solver.hpp>

namespace itl { namespace pc {

/// Class for applying \tparam Preconditioner only on a sub-matrix
/** Other entries are just copied. 
    Optionally preconditioner can be referred from outside instead of storing it by setting
    \tparam Store to true. 
**/
template <typename Preconditioner, typename Matrix, bool Store= true>
class sub_matrix_pc
{
    typedef mtl::dense_vector<bool>                                               tag_type;
    typedef typename boost::mpl::if_c<Store, Preconditioner, const Preconditioner&>::type pc_type;

    struct matrix_container
    {
	matrix_container() : Ap(0) {} 

	matrix_container(const tag_type& tags, const Matrix& src)
	{
	    using std::size_t; 
	    using namespace mtl;

	    mtl::dense_vector<size_t> perm(size(tags));
	    size_t n= 0;
	    for (size_t i= 0; i < size(tags); ++i) {
		perm[i]= n;
		if (tags[i])
		    ++n;
	    }

	    Ap= new Matrix(n, n);

	    typename traits::row<Matrix>::type             row(src); 
	    typename traits::col<Matrix>::type             col(src); 
	    typename traits::const_value<Matrix>::type     value(src); 
	    typedef typename traits::range_generator<tag::major, Matrix>::type  cursor_type;

	    mat::inserter<Matrix> ins(*Ap, Ap->nnz() / Ap->dim1());
	    
	    for (cursor_type cursor = mtl::begin<tag::major>(src), cend = mtl::end<tag::major>(src); 
		 cursor != cend; ++cursor) {
		// std::cout << dest << '\n';
	    
		typedef typename traits::range_generator<tag::nz, cursor_type>::type icursor_type;
		for (icursor_type icursor = mtl::begin<tag::nz>(cursor), icend = mtl::end<tag::nz>(cursor); 
		     icursor != icend; ++icursor) 
		    if (tags[row(*icursor)] && tags[col(*icursor)])
			ins(perm[row(*icursor)], perm[col(*icursor)]) << value(*icursor); 
	    }
	}

	~matrix_container() { delete Ap; }

	Matrix* Ap;
    };

    size_t count_entries() const
    {
	using mtl::size;
	size_t n= 0;
	for (size_t i= 0; i < size(tags); ++i) 
	    if (tags[i]) ++n;
	return n;
    }

  public:

    sub_matrix_pc(const tag_type& tags, const Matrix& A)
      : tags(tags), n(count_entries()), mc(tags, A), P(*mc.Ap)
    {
	BOOST_STATIC_ASSERT((Store));
	delete mc.Ap; 
	mc.Ap= 0;
    }

#if 0 // to do later
    sub_matrix_pc(const tag_type& tags, const Preconditioner& P)
      : tags(tags), P(P)
    {
	// check sizes
    }
#endif

  private:
    template <typename VectorIn>
    VectorIn& create_x0(VectorIn) const
    {
	static VectorIn  x0(n);
	return x0;
    }

    template <typename VectorOut>
    VectorOut& create_y0(VectorOut) const
    {
	static VectorOut  y0(n);
	return y0;
    }

    template <typename VectorIn>
    void restrict(const VectorIn& x, VectorIn& x0) const
    {
	for (size_t i= 0, j= 0; i < size(tags); ++i) 
	    if (tags[i]) 
		x0[j++]= x[i];
    }

    template <typename VectorIn, typename VectorOut>
    void prolongate(const VectorIn& x, const VectorOut& y0, VectorOut& y) const
    {
	for (size_t i= 0, j= 0; i < size(tags); ++i) 
	    if (tags[i]) 
		y[i]= y0[j++];
	    else
		y[i]= x[i];	
    }

  public:
    /// Solve Px = y approximately on according sub-system; remaining entries are copied
    template <typename VectorIn, typename VectorOut>
    void solve(const VectorIn& x, VectorOut& y) const
    {
	mtl::vampir_trace<5056> tracer;
	y.checked_change_resource(x);

	VectorIn&  x0= create_x0(x);
	VectorOut& y0= create_y0(y);

	restrict(x, x0);
	P.solve(x0, y0);
	// y0= solve(P, x0); // doesn't compile yet for unknown reasons
	prolongate(x, y0, y);
    }

    /// Solve Px = y approximately on according sub-system; remaining entries are copied
    template <typename VectorIn, typename VectorOut>
    void adjoint_solve(const VectorIn& x, VectorOut& y) const
    {
	mtl::vampir_trace<5057> tracer;
	y.checked_change_resource(x);

	VectorIn&  x0= create_x0(x);
	VectorOut& y0= create_y0(y);

	restrict(x, x0);
	P.adjoint_solve(x0, y0);
	// y0= adjoint_solve(P, x0); // doesn't compile yet for unknown reasons
	prolongate(x, y0, y);
    }

  private:
    tag_type          tags;
    size_t            n;
    matrix_container  mc;
    pc_type           P;
};

template <typename Preconditioner, typename Matrix, bool Store, typename Vector>
solver<sub_matrix_pc<Preconditioner, Matrix, Store>, Vector, false>
inline solve(const sub_matrix_pc<Preconditioner, Matrix, Store>& P, const Vector& x)
{
    return solver<sub_matrix_pc<Preconditioner, Matrix, Store>, Vector, false>(P, x);
}

template <typename Preconditioner, typename Matrix, bool Store, typename Vector>
solver<sub_matrix_pc<Preconditioner, Matrix, Store>, Vector, true>
inline adjoint_solve(const sub_matrix_pc<Preconditioner, Matrix, Store>& P, const Vector& x)
{
    return solver<sub_matrix_pc<Preconditioner, Matrix, Store>, Vector, true>(P, x);
}


}} // namespace itl::pc

#endif // ITL_PC_SUB_MATRIX_PC_INCLUDE
