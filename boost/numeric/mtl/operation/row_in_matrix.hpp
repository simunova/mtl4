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

#ifndef MTL_ROW_IN_MATRIX_INCLUDE
#define MTL_ROW_IN_MATRIX_INCLUDE

#include <boost/mpl/if.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

/// Type of row in matrix as vector if accessible
template <typename Matrix>
struct RowInMatrix 
{
    static const bool exists= false;
};

template <typename Value, typename Parameters>
struct RowInMatrix<mtl::mat::dense2D<Value, Parameters> > 
{
    typedef mtl::mat::dense2D<Value, Parameters> ref_type;
    typedef typename ref_type::size_type       size_type;
    typedef typename ref_type::value_type      value_type;
    typedef mtl::vec::parameters<row_major>    vec_para;

    static const bool aligned= boost::is_same<typename Parameters::orientation, row_major>::value;
    static const bool exists= true;

    typedef typename boost::mpl::if_c<
	aligned
      , vec::dense_vector<Value, vec_para> 
      , vec::strided_vector_ref<Value, vec_para> 
    >::type type;

    static inline type apply(ref_type& A, size_type row, const irange& col_range)
    {
	return dispatch<type, ref_type>(A, row, col_range, boost::mpl::bool_<aligned>());
    }

    template <typename Ref>
    static inline size_type vector_size(const Ref& A, const irange& col_range)
    {
	using std::min;
	
	vampir_trace<1003> tracer;
	size_type finish= min(col_range.finish(), num_cols(A));
	return col_range.start() < finish ? finish - col_range.start() : 0;
    }

    template <typename Return, typename Ref>
    static inline Return dispatch(Ref& A, size_type row, const irange& col_range, boost::mpl::true_)
    {
	vampir_trace<2023> tracer;
	return Return(vector_size(A, col_range), const_cast<value_type*>(&A[row][col_range.start()])); // TODO make work without const cast
    }

    template <typename Return, typename Ref>
    static inline Return dispatch(Ref& A, size_type row, const irange& col_range, boost::mpl::false_)
    {
	vampir_trace<1004> tracer;
	return Return(vector_size(A, col_range), &A[row][col_range.start()], num_rows(A));
    }	 
};

template <typename Value, typename Parameters>
struct RowInMatrix<const mtl::mat::dense2D<Value, Parameters> > 
{
    typedef mtl::mat::dense2D<Value, Parameters> const   ref_type;
    typedef mtl::mat::dense2D<Value, Parameters>         ref2_type;
    typedef typename ref2_type::size_type      size_type;
    typedef vec::parameters<row_major>      vec_para;

    static const bool aligned= boost::is_same<typename Parameters::orientation, row_major>::value;
    static const bool exists= true;

    typedef typename boost::mpl::if_c<
        aligned
      , vec::dense_vector<Value, vec_para> // TODO needs constification !!!
      , vec::strided_vector_ref<const Value, vec_para> 
    >::type type;

    static inline type apply(ref_type& A, size_type row, const irange& col_range)
    {
	return RowInMatrix<ref2_type>::template dispatch<type, ref_type>(A, row, col_range, boost::mpl::bool_<aligned>());
    }
};


} // namespace mtl

#endif // MTL_ROW_IN_MATRIX_INCLUDE
