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

#ifndef ITL_PC_ILU_0_INCLUDE
#define ITL_PC_ILU_0_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/invert_diagonal.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/matrix/strict_lower.hpp>
#include <boost/numeric/mtl/matrix/upper.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include <boost/numeric/itl/pc/ilu.hpp>

namespace itl { namespace pc {

// Dummy type to perform factorization in initializer list not in 
struct ilu_0_factorizer
{
    template <typename Matrix, typename L_type, typename U_type>
    ilu_0_factorizer(const Matrix &A, L_type& L, U_type& U)
    {   factorize(A, L, U, mtl::traits::is_sparse<Matrix>());  }

    template <typename Matrix, typename L_type, typename U_type>
    void factorize(const Matrix&, L_type&, U_type&, boost::mpl::false_)
    {  MTL_THROW_IF(true, mtl::logic_error("ILU is not intended for dense matrices")); }

    template <typename Matrix, typename L_type, typename U_type>
    void factorize(const Matrix& A, L_type& L, U_type& U, boost::mpl::true_)
    {
	using namespace mtl; using namespace mtl::tag;  using mtl::traits::range_generator;  
	using math::reciprocal; 
	MTL_THROW_IF(num_rows(A) != num_cols(A), mtl::matrix_not_square());
	mtl::vampir_trace<5038> tracer;

	typedef typename mtl::Collection<Matrix>::value_type      value_type;
	typedef typename mtl::Collection<Matrix>::size_type       size_type;
	typedef mtl::mat::parameters<mtl::row_major, mtl::index::c_index, mtl::non_fixed::dimensions, false, size_type> para;
	typedef mtl::mat::compressed2D<value_type, para>  LU_type;
	LU_type LU(A);

	typedef typename range_generator<row, LU_type>::type      cur_type;    
	typedef typename range_generator<nz, cur_type>::type      icur_type;            
	typename mtl::traits::col<LU_type>::type                  col(LU);
	typename mtl::traits::value<LU_type>::type                value(LU); 
	mtl::vec::dense_vector<value_type, mtl::vec::parameters<> > inv_dia(num_rows(A));
	cur_type ic= begin<row>(LU), iend= end<row>(LU);
	for (size_type i= 0; ic != iend; ++ic, ++i) {

	    for (icur_type kc= begin<nz>(ic), kend= end<nz>(ic); kc != kend; ++kc) {
		size_type k= col(*kc);
		if (k == i) break;

		value_type aik= value(*kc) * inv_dia[k];
		value(*kc, aik);

		for (icur_type jc= kc + 1; jc != kend; ++jc)
		    value(*jc, value(*jc) - aik * LU[k][col(*jc)]);
		// std::cout << "LU after eliminating A[" << i << "][" << k << "] =\n" << LU;			  
	    }
	    inv_dia[i]= reciprocal(LU[i][i]);
	}
	invert_diagonal(LU); 
	L= strict_lower(LU);
	U= upper(LU);
    }  
};

template <typename Matrix, typename Value= typename mtl::Collection<Matrix>::value_type>
class ilu_0
  : public ilu<Matrix, ilu_0_factorizer, Value>
{
    typedef ilu<Matrix, ilu_0_factorizer, Value> base;
  public:
    ilu_0(const Matrix& A) : base(A) {}
};

// ic_0_evaluator not needed IC(0) and ILU(0) do the same at the upper triangle ;-)

}} // namespace itl::pc


#endif // ITL_PC_ILU_0_INCLUDE
