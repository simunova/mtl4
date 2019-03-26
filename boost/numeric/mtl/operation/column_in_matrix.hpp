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

#ifndef MTL_COLUMN_IN_MATRIX_INCLUDE
#define MTL_COLUMN_IN_MATRIX_INCLUDE

#include <boost/mpl/if.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>

namespace mtl {


/// Type of column in matrix as vector if accessible
template <typename Matrix>
struct ColumnInMatrix
{
    static const bool exists= false;
};

template <typename Value, typename Parameters>
struct ColumnInMatrix<mtl::mat::dense2D<Value, Parameters> > 
{
    typedef mtl::mat::dense2D<Value, Parameters>         ref_type;
    typedef typename ref_type::size_type       size_type;
    typedef typename ref_type::value_type      value_type;

    static const bool aligned= !boost::is_same<typename Parameters::orientation, row_major>::value;
    static const bool exists= true;

    typedef typename boost::mpl::if_c<
	aligned
      , vec::dense_vector<Value, vec::parameters<> > 
      , vec::strided_vector_ref<Value, vec::parameters<> > 
    >::type type;

    static inline type apply(ref_type& A, const irange& row_range, size_type col)
    {
	return dispatch<type, ref_type>(A, row_range, col, boost::mpl::bool_<aligned>());
    }

    template <typename Ref>
    static inline size_type vector_size(const Ref& A, const irange& row_range)
    {
	using std::min;
	size_type finish= min(row_range.finish(), num_rows(A));
	return row_range.start() < finish ? finish - row_range.start() : 0;
    }

    template <typename Return, typename Ref>
    static inline Return dispatch(Ref& A, const irange& row_range, size_type col, boost::mpl::true_)
    {
	return Return(vector_size(A, row_range), const_cast<value_type*>(&A[row_range.start()][col])); // TODO make work without const cast
    }

    template <typename Return, typename Ref>
    static inline Return dispatch(Ref& A, const irange& row_range, size_type col, boost::mpl::false_)
    {
	return Return(vector_size(A, row_range), &A[row_range.start()][col], num_cols(A));
    }	 
};

template <typename Value, typename Parameters>
struct ColumnInMatrix<const mtl::mat::dense2D<Value, Parameters> > 
{
    typedef mtl::mat::dense2D<Value, Parameters> const   ref_type;
    typedef mtl::mat::dense2D<Value, Parameters>         ref2_type;
    typedef typename ref2_type::size_type      size_type;

    static const bool aligned= !boost::is_same<typename Parameters::orientation, row_major>::value;
    static const bool exists= true;

    typedef typename boost::mpl::if_c<
	aligned
      , vec::dense_vector<Value, vec::parameters<> > // TODO needs constification !!!
      , vec::strided_vector_ref<const Value, vec::parameters<> > 
    >::type type;

    static inline type apply(ref_type& A, const irange& row_range, size_type col)
    {
	return ColumnInMatrix<ref2_type>::template dispatch<type, ref_type>(A, row_range, col, boost::mpl::bool_<aligned>());
    }
};
    
} // namespace mtl

#endif // MTL_COLUMN_IN_MATRIX_INCLUDE
