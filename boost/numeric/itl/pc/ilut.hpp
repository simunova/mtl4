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

#ifndef ITL_PC_ILUT_INCLUDE
#define ITL_PC_ILUT_INCLUDE

#include <boost/numeric/mtl/vector/sparse_vector.hpp>
#include <boost/numeric/mtl/operation/invert_diagonal.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace itl { namespace pc {


struct ilut_factorizer
{

    template <typename Matrix, typename Para, typename L_type, typename U_type>
    ilut_factorizer(const Matrix &A, const Para& p, L_type& L, U_type& U)
    { factorize(A, p, L, U, mtl::traits::is_row_major<Matrix>()); }

    // column-major matrices are copied first
    template <typename Matrix, typename Para, typename L_type, typename U_type, bool B>
    void factorize(const Matrix &A, const Para& p, L_type& L, U_type& U, boost::mpl::bool_<B>)
    {
	typedef typename mtl::Collection<Matrix>::value_type      value_type;
	typedef typename mtl::Collection<Matrix>::size_type       size_type;
	typedef mtl::mat::parameters<mtl::row_major, mtl::index::c_index, mtl::non_fixed::dimensions, false, size_type> para;
	typedef mtl::mat::compressed2D<value_type, para>  LU_type;
	LU_type LU(A);
	factorize(LU, p, L, U, boost::mpl::true_());
    }

    // According Yousef Saad: ILUT, NLAA, Vol 1(4), 387-402 (1994)
#if 0
    template <typename Value, typename MPara, typename Para, typename L_type, typename U_type>
    factorize(const mtl::mat::compressed2D<Value, MPara>& A, const Para& p, L_type& L, U_type& U, boost::mpl::true_)
#endif

    template <typename Matrix, typename Para, typename L_type, typename U_type>
    void factorize(const Matrix& A, const Para& p, L_type& L, U_type& U, boost::mpl::true_)

    {   
	mtl::vampir_trace<5049> tracer;
	using std::abs; using mtl::traits::range_generator; using mtl::begin; using mtl::end;
	using namespace mtl::tag;
	MTL_THROW_IF(num_rows(A) != num_cols(A), mtl::matrix_not_square());

	typedef typename mtl::Collection<Matrix>::value_type      value_type;
	typedef typename mtl::Collection<Matrix>::size_type       size_type;
	typedef typename range_generator<row, Matrix>::type       cur_type;    
	typedef typename range_generator<nz, cur_type>::type      icur_type;            
	typename mtl::traits::col<Matrix>::type                   col(A);
	typename mtl::traits::const_value<Matrix>::type           value(A); 

	size_type n= num_rows(A);
	L.change_dim(n, n); 
	U.change_dim(n, n);
	{
	    mtl::mat::inserter<L_type> L_ins(L, p.first);
	    mtl::mat::inserter<U_type> U_ins(U, p.first + 1); // plus one for diagonal
	
	    mtl::sparse_vector<value_type> vec(n); // corr. row in paper
	    cur_type ic= begin<row>(A); // , iend= end<row>(A);
	    for (size_type i= 0; i < n; ++i, ++ic) {
		
		for (icur_type kc= begin<nz>(ic), kend= end<nz>(ic); kc != kend; ++kc) // row= A[i][*]
		    vec.insert(col(*kc), value(*kc));
		// std::cerr << "vec_" << i << " = " << vec << std::endl;
		value_type tau_i= p.second * two_norm(vec); // threshold for i-th row
		// loop over non-zeros in vec; changes in vec considered
		for (size_type j= 0; j < vec.nnz() && vec.index(j) < i; j++) {
		    size_type k= vec.index(j);
		    value_type ukk= U_ins.value(k, k);
		    MTL_DEBUG_THROW_IF(ukk == value_type(0), mtl::missing_diagonal());
		    value_type vec_k= vec.value(j)/= ukk;
		    // std::cout << "vec after updating from U[" << k << "][" << k << "] is " << vec << '\n';
		    for (size_type j0= U_ins.ref_major()[k], j1= U_ins.ref_slot_ends()[k]; j0 < j1; j0++) { // U[k][k+1:n]
			size_type k1= U_ins.ref_minor()[j0];
			if (k1 > k)
			    vec[k1]-= vec_k * U_ins.ref_elements()[j0];
			// std::cout << "vec after updating from U[" << k << "][" << k1 << "] is " << vec << '\n';
		    }
		    // if (i > 1000 && i < 1010) std::cout << "vec before crop in row " << i << ", updating from row " << k << ": \n" << vec << "\n";
		    vec.crop(tau_i);
		    // if (i > 1000 && i < 1010) std::cout << "vec after crop: \n" << vec << "\n";
		}
		// std::cerr << "vec_" << i << " = " << vec << std::endl;
		vec.sort_on_data();
		// std::cerr << "vec_" << i << " sorted on data = " << vec << std::endl;
		
		// std::cout << "vec at " << i << " is " << vec << '\n';
		// mtl::vampir_trace<9904> tracer2;
		bool diag_found= false;
		for (size_type cntu= 0, cntl= 0, j= 0; j < vec.nnz() && (cntu < p.first || cntl < p.first); j++) {
		    size_type k= vec.index(j);
		    value_type v= vec.value(j);
		    // if (abs(v) < tau_i) break;
		    if (i == k) {
			U_ins[i][i] << v; diag_found= true;
		    } else if (i < k) {
			if (cntu++ < p.first)
			    U_ins[i][k] << v;
		    } else // i > k
			if (cntl++ < p.first)
			    L_ins[i][k] << v;		
		}
		if (!diag_found) std::cerr << "Deleted diagonal!!!!\n";
		vec.make_empty();
	    }
	} // destroy inserters
	invert_diagonal(U);
    }
};

// Not usable yet !!!!!
template <typename Matrix, typename Value= typename mtl::Collection<Matrix>::value_type>
class ilut
  : public ilu<Matrix, ilut_factorizer, Value>
{
    typedef ilu<Matrix, ilut_factorizer, Value> base;
  public:
    ilut(const Matrix& A, std::size_t p, typename mtl::Collection<Matrix>::value_type tau) 
      : base(A, std::make_pair(p, tau)) {}
};

}} // namespace itl::pc

#endif // ITL_PC_ILUT_INCLUDE
