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

#ifndef MTL_CRTP_BASE_MATRIX_INCLUDE
#define MTL_CRTP_BASE_MATRIX_INCLUDE

#include <iostream>
#include <algorithm>
#include <boost/mpl/bool.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/numeric/mtl/operation/print.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/matrix_bracket.hpp>
#include <boost/numeric/mtl/operation/copy.hpp>
#include <boost/numeric/mtl/operation/mult.hpp>
#include <boost/numeric/mtl/operation/right_scale_inplace.hpp>
#include <boost/numeric/mtl/operation/divide_by_inplace.hpp>
#include <boost/numeric/mtl/matrix/all_mat_expr.hpp>
#include <boost/numeric/mtl/matrix/diagonal_setup.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/assert.hpp>
#include <boost/numeric/mtl/utility/eval_dense.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/utility/iset.hpp>
#include <boost/numeric/mtl/operation/mult_assign_mode.hpp>
#include <boost/numeric/mtl/operation/compute_factors.hpp>
#include <boost/numeric/mtl/operation/column_in_matrix.hpp>
#include <boost/numeric/mtl/operation/row_in_matrix.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#ifdef MTL_WITH_INITLIST
# include <initializer_list>
#endif

namespace mtl { namespace mat {

template <typename Source, typename Matrix>
struct crtp_assign 
{
    Matrix& operator()(const Source& source, Matrix& matrix)
    {
	return assign(source, matrix, typename ashape::ashape<Source>::type());
    }
private:
    /// Assign scalar to a matrix by setting the matrix to a multiple of unity matrix
    /** Uses internally \sa diagonal_setup, for details see there. **/
    Matrix& assign(const Source& source, Matrix& matrix, ashape::scal)
    {
	vampir_trace<3055> tracer;
	MTL_CRASH_IF(num_rows(matrix) * num_cols(matrix) == 0, 
		     "Trying to initialize a 0 by 0 matrix with a value");
	diagonal_setup(matrix, source);
	return matrix;
    }

    /// Assign matrix expressions by copying except for some special expressions
    Matrix& assign(const Source& source, Matrix& matrix, typename ashape::ashape<Matrix>::type)
    {
	vampir_trace<3056> tracer;
	// Self-assignment between different types shouldn't happen.	
	matrix.checked_change_resource(source);
	matrix_copy(source, matrix);
	return matrix;
    }
};



/// Assign sum by assigning first argument and adding second
/*  Note that this is more special then assigning arbitrary expressions including matrices itself
    because mat_mat_plus_expr <E1, E2> is a derived class from mat_expr < MatrixSrc >. **/
template <typename E1, typename E2, typename Matrix>
struct crtp_assign<mat_mat_plus_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_plus_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3056> tracer;
	matrix.checked_change_resource(src.first);
	matrix= src.first;
	return matrix+= src.second;
    }
};

/// Assign difference by assigning first argument and subtracting second
/*  Note that this is more special then assigning arbitrary expressions including matrices itself
    because mat_mat_minus_expr <E1, E2> is a derived class from mat_expr < MatrixSrc >. **/
template <typename E1, typename E2, typename Matrix>
struct crtp_assign<mat_mat_minus_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_minus_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3057> tracer;
	matrix.checked_change_resource(src.first);
	matrix= src.first;
	return matrix-= src.second;
    }
};

/// Assign product by calling mult
template <typename E1, typename E2, typename Matrix>
struct crtp_assign<mat_mat_times_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_times_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<4012> tracer;
	operation::compute_factors<Matrix, mat_mat_times_expr<E1, E2> > factors(src);
	matrix.checked_change_resource(factors.first, factors.second);
	mult(factors.first, factors.second, matrix);
	return matrix;
    }
}; 


/// Assign element-wise product 
template <typename E1, typename E2, typename Matrix>
struct crtp_assign<mat_mat_ele_times_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_ele_times_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3028> tracer;
	operation::compute_factors<Matrix, mat_mat_ele_times_expr<E1, E2> > factors(src);
	matrix.checked_change_resource(factors.first);
	matrix= factors.first;
	return matrix.ele_rscale(factors.second);
    }
}; 



/// Assign c-style 2D-array, because it's easier to initialize.
template <typename Value, unsigned Rows, unsigned Cols, typename Matrix>
struct crtp_assign<Value[Rows][Cols], Matrix>
{
    Matrix& operator()(const Value src[Rows][Cols], Matrix& matrix)
    {
	vampir_trace<3059> tracer;
	typedef typename Collection<Matrix>::size_type size_type;

	matrix.checked_change_dim(Rows, Cols);
	inserter<Matrix>  ins(matrix, matrix.dim2());
	
	for (size_type r= 0; r < Rows; ++r)
	    for (size_type c= 0; c < Cols; ++c)
		ins(r, c) << src[r][c];
	return matrix;
    }
};

#if defined(MTL_WITH_INITLIST) && defined(MTL_WITH_AUTO) && defined(MTL_WITH_RANGEDFOR)
    /// Constructor for initializer list \p values 
    template <typename Value2, typename Matrix>
    struct crtp_assign<std::initializer_list<std::initializer_list<Value2> >, Matrix>
    {
	Matrix& operator()(std::initializer_list<std::initializer_list<Value2> > values, Matrix& matrix)
	{
	    typedef typename Collection<Matrix>::size_type size_type;
	    size_type nr= values.size(), nc= nr > 0? values.begin()->size() : 0;
	    matrix.checked_change_dim(nr, nc);
	    inserter<Matrix>  ins(matrix, matrix.dim2());

	    size_t r= 0;
	    for (auto l : values) {
		size_t c= 0;	    
		MTL_CRASH_IF(l.size() != nc, "All sub-lists must have same size!");
		for (auto v : l)
		    ins(r, c++) << v;
		r++;
	    }
	    return matrix;
	}
    };
#endif


template <typename Vector, typename Matrix>
struct crtp_assign<multi_vector<Vector>, Matrix>
{
    Matrix& operator()(const multi_vector<Vector>& src, Matrix& matrix)
    {
	vampir_trace<3060> tracer;
	typedef typename Collection<Matrix>::size_type size_type;

	matrix.checked_change_resource(src);
	// del checked_change_dim(num_rows(src), num_cols(src));
	inserter<Matrix>  ins(matrix);
	
	for (size_type r= 0; r < num_rows(src); ++r)
	    for (size_type c= 0; c < num_cols(src); ++c)
		ins(r, c) << src[r][c];
	return matrix;
    }
};

template <typename Matrix>
struct crtp_assign<identity2D, Matrix>
{
    Matrix& operator()(const identity2D& src, Matrix& matrix)
    {
	typedef typename Collection<Matrix>::size_type  size_type;
	typedef typename Collection<Matrix>::value_type value_type;

	matrix.checked_change_resource(src);
	{
	    inserter<Matrix>  ins(matrix);	
	    for (size_type r= 0; r < num_rows(src); ++r)
		for (size_type c= 0; c < num_cols(src); ++c)
		    ins(r, c) << value_type(int(r == c));
	}
	return matrix;
    }
};




/// Assign content of a file to the matrix
template <typename IFStream, typename OFStream, typename Matrix>
struct crtp_assign<io::matrix_file<IFStream, OFStream>, Matrix>
{
    Matrix& operator()(const io::matrix_file<IFStream, OFStream>& file, Matrix& matrix)
    {
	vampir_trace<3029> tracer;
	IFStream stream(file.file_name().c_str());
	stream >> matrix;
	return matrix;
    }
};
	
/// Assign-add matrix expressions by incrementally copying except for some special expressions
template <typename Source, typename Matrix>
struct crtp_plus_assign 
{
    Matrix& operator()(const Source& source, Matrix& matrix)
    {
	vampir_trace<3030> tracer;
	return assign(source, matrix, typename ashape::ashape<Source>::type());
    }
  private:
    Matrix& assign(const Source& source, Matrix& matrix, typename ashape::ashape<Matrix>::type)
    {
	matrix_copy_plus(source, matrix);
	return matrix;
    }
};

/// Assign-add sum by adding both arguments
/** Note that this is more special then assigning arbitrary expressions including matrices itself
	because mat_mat_plus_expr <E1, E2> is a derived class from 
	mat_expr < MatrixSrc >. **/
template <typename E1, typename E2, typename Matrix>
struct crtp_plus_assign<mat_mat_plus_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_plus_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3030> tracer;
	matrix+= src.first;
	return matrix+= src.second;
    }
};

template <typename E1, typename E2, typename Matrix>
struct crtp_plus_assign<mat_mat_minus_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_minus_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3030> tracer;
	matrix+= src.first;
	return matrix-= src.second;
    }
};

template <typename E1, typename E2, typename Matrix>
struct crtp_plus_assign<mat_mat_ele_times_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_ele_times_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3030> tracer;
	Matrix Prod(ele_prod(src.first, src.second));
	return matrix+= Prod;
    }
};

template <typename E1, typename E2, typename Matrix>
struct crtp_plus_assign<mat_mat_times_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_times_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3030> tracer;
	operation::compute_factors<Matrix, mat_mat_times_expr<E1, E2> > factors(src);
	gen_mult(factors.first, factors.second, matrix, assign::plus_sum(), 
		 tag::flat<tag::matrix>(), tag::flat<tag::matrix>(), tag::flat<tag::matrix>());
	return matrix;
    }
};


/// Assign-subtract matrix expressions by decrementally copying except for some special expressions
template <typename Source, typename Matrix>
struct crtp_minus_assign 
{
    Matrix& operator()(const Source& source, Matrix& matrix)
    {
	vampir_trace<3031> tracer;
	return assign(source, matrix, typename ashape::ashape<Source>::type());
    }
private:
    Matrix& assign(const Source& source, Matrix& matrix, typename ashape::ashape<Matrix>::type)
    {
	matrix_copy_minus(source, matrix);
	return matrix;
    }
};

/// Assign-subtract sum by adding both arguments
/** Note that this is more special then assigning arbitrary expressions including matrices itself
	because mat_mat_plus_expr <E1, E2> is a derived class from 
	mat_expr < MatrixSrc >. **/
template <typename E1, typename E2, typename Matrix>
struct crtp_minus_assign<mat_mat_plus_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_plus_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3031> tracer;
	matrix-= src.first;
	return matrix-= src.second;
    }
};

/// Assign-subtracting difference by subtracting first argument and adding the second one
/** Note that this is more special then assigning arbitrary expressions including matrices itself
	because mat_mat_minus_expr <E1, E2> is a derived class from 
	mat_expr < MatrixSrc >. **/
template <typename E1, typename E2, typename Matrix>
struct crtp_minus_assign<mat_mat_minus_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_minus_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3031> tracer;
	matrix-= src.first;
	return matrix+= src.second;
    }
};

template <typename E1, typename E2, typename Matrix>
struct crtp_minus_assign<mat_mat_ele_times_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_ele_times_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3031> tracer;
	Matrix Prod(ele_prod(src.first, src.second));
	return matrix-= Prod;
    }
};

/// Assign-subtract product by calling gen_mult
/** Note that this does not work for arbitrary expressions. **/
template <typename E1, typename E2, typename Matrix>
struct crtp_minus_assign<mat_mat_times_expr<E1, E2>, Matrix> 
{
    Matrix& operator()(const mat_mat_times_expr<E1, E2>& src, Matrix& matrix)
    {
	vampir_trace<3031> tracer;
	operation::compute_factors<Matrix, mat_mat_times_expr<E1, E2> > factors(src);
	gen_mult(factors.first, factors.second, matrix, assign::minus_sum(), tag::flat<tag::matrix>(), tag::flat<tag::matrix>(), tag::flat<tag::matrix>());
	return matrix;
    }
};



/// Base class to provide matrix assignment operators generically 
template <typename Matrix, typename ValueType, typename SizeType>
struct crtp_matrix_assign
{
private:

    // For (compatible) dense matrices do a loop over all entries
    template <typename Source>
    Matrix& density_assign(const Source& src, boost::mpl::true_)
    {
	vampir_trace<3032> tracer;
	// typedef typename Collection<Source>::size_type size_type;
	typedef unsigned size_type;

	// std::cout << "Dense assignment\n";
	checked_change_resource(src);

	Matrix& matrix= static_cast<Matrix&>(*this);
	for (size_type r= 0; r < num_rows(matrix); ++r)
	    for (size_type c= 0; c < num_cols(matrix); ++c)
		matrix[r][c]= src[r][c];
	return matrix;
    }

    // If sparse matrices are involved evaluate step-wise (or assignment from scalar)
    template <typename Source>
    Matrix& density_assign(const Source& src, boost::mpl::false_)
    {
	// std::cout << "Sparse assignment\n";
	return crtp_assign<Source, Matrix>()(src, static_cast<Matrix&>(*this));
    }

    
    // For (compatible) dense matrices do a loop over all entries
    template <typename Source>
    Matrix& density_plus_assign(const Source& src, boost::mpl::true_)
    {
	vampir_trace<3033> tracer;
	// typedef typename Collection<Source>::size_type size_type;
	typedef unsigned size_type;

	// std::cout << "Dense assignment\n";
	checked_change_resource(src);
	// del checked_change_dim(num_rows(src), num_cols(src));

	Matrix& matrix= static_cast<Matrix&>(*this);
	for (size_type r= 0; r < num_rows(matrix); ++r)
	    for (size_type c= 0; c < num_cols(matrix); ++c)
		matrix[r][c]+= src[r][c];
	return matrix;
    }

    // If sparse matrices are involved evaluate step-wise (or assignment from scalar)
    template <typename Source>
    Matrix& density_plus_assign(const Source& src, boost::mpl::false_)
    {
	// std::cout << "Sparse assignment\n";
	return crtp_plus_assign<Source, Matrix>()(src, static_cast<Matrix&>(*this));
    }

    // For (compatible) dense matrices do a loop over all entries
    template <typename Source>
    Matrix& density_minus_assign(const Source& src, boost::mpl::true_)
    {
	vampir_trace<3034> tracer;
	// typedef typename Collection<Source>::size_type size_type;
	typedef unsigned size_type;

	// std::cout << "Dense assignment\n";
	checked_change_resource(src);

	Matrix& matrix= static_cast<Matrix&>(*this);
	for (size_type r= 0; r < num_rows(matrix); ++r)
	    for (size_type c= 0; c < num_cols(matrix); ++c)
		matrix[r][c]-= src[r][c];
	return matrix;
    }

    // If sparse matrices are involved evaluate step-wise (or assignment from scalar)
    template <typename Source>
    Matrix& density_minus_assign(const Source& src, boost::mpl::false_)
    {
	// std::cout << "Sparse assignment\n";
	return crtp_minus_assign<Source, Matrix>()(src, static_cast<Matrix&>(*this));
    }

    // For (compatible) dense matrices do a loop over all entries
    template <typename Source>
    Matrix& density_ele_rscale(const Source& src, boost::mpl::true_)
    {
	vampir_trace<3035> tracer;
	// typedef typename Collection<Source>::size_type size_type;
	typedef unsigned size_type;

	// std::cout << "Dense assignment\n";
	checked_change_resource(src);
	// del checked_change_dim(num_rows(src), num_cols(src));

	Matrix& matrix= static_cast<Matrix&>(*this);
	for (size_type r= 0; r < num_rows(matrix); ++r)
	    for (size_type c= 0; c < num_cols(matrix); ++c)
		matrix[r][c]*= src[r][c];
	return matrix;
    }

    // If sparse matrices are involved evaluate step-wise (or assignment from scalar)
    template <typename Factor>
    Matrix& density_ele_rscale(const Factor& alpha, boost::mpl::false_)
    {
	// std::cout << "Sparse assignment\n";
	matrix_copy_ele_times(alpha, static_cast<Matrix&>(*this));
	return static_cast<Matrix&>(*this);
    }

  public:

    /// Check wether source and target have compatible resources, generalization of check_dim
    /** For expressions like A= B + C, A can be set to the size of B and C if still is 0 by 0. **/
    template <typename Src>
    void check_resource(const Src& src) const 
    {	check_resource(src, typename mtl::traits::category<Matrix>::type());    }

    // Default case just check_dim
    template <typename Src>
    void check_resource(const Src& src, tag::universe) const 
    {	check_dim(num_rows(src), num_cols(src));    }

    /// Check wether source and target have compatible resources and wether target has already resources
    /** For expressions like A+= B + C, A must be already larger then 0 by 0 and compatible to B and C. **/
    //  Generalization with 2 arguments might be needed (check rows from first and columns from second)
    template <typename Src>
    void check_ready_resource(const Src& src) const 
    {
	MTL_CRASH_IF(num_rows(src) * num_cols(src) == 0, "Need non-empty matrix!");
	check_resource(src);
    }

    /// Check wether source and target have compatible resources and adapt empty target
    /** For expressions like A= B + C, A can be set to the size of B and C if still is 0 by 0. **/
    template <typename Src>
    void checked_change_resource(const Src& src) 
    {	checked_change_resource(src, src);   }

    /// Check whether source and target have compatible resources and adapt empty target
    /** For expressions like A= B + C, A can be set to the size of B and C if still is 0 by 0. **/
    template <typename Src1, typename Src2>
    void checked_change_resource(const Src1& src1, const Src2& src2)
    {   checked_change_resource_aux(src1, src2, typename mtl::traits::category<Matrix>::type());    }

    template <typename Src1, typename Src2>
    void checked_change_resource_aux(const Src1& src1, const Src2& src2, tag::universe) 
    {   checked_change_dim(SizeType(num_rows(src1)), SizeType(num_cols(src2)));  }


    /// Check whether matrix sizes are compatible or if matrix is 0 by 0 change it to r by c.
    /** Deprecated, superseded by checked_change_resource. **/ 
    void checked_change_dim(SizeType r, SizeType c)
    {
	Matrix& matrix= static_cast<Matrix&>(*this);
	matrix.check_dim(r, c);
	matrix.change_dim(r, c);
    }

    /// Templated assignment implemented by functor to allow for partial specialization
    // Despite there is only an untemplated assignment and despite the disable_if MSVC whines about ambiguity :-!
    // Scalar assignment is also taking out because it has another return type
    template <typename Source>
    typename boost::disable_if_c<boost::is_same<Matrix, Source>::value 
                                   || boost::is_same<typename ashape::ashape<Source>::type, ashape::scal>::value,
				 Matrix&>::type
    operator=(const Source& src)
    {
	return density_assign(src, boost::mpl::bool_< boost::is_same<typename ashape::ashape<Matrix>::type, 
			                                             typename ashape::ashape<Source>::type>::value 
			                              && mtl::traits::eval_dense< mat_mat_asgn_expr<Matrix, Source> >::value >());
    }

    // Helper type for assigning scalars to handle both A= a; and A= a, b, c;
    template <typename Source>
    struct scalar_assign 
    {
	scalar_assign(Source src, Matrix& matrix) 
	  : src(src), with_comma(false), r(0), c(0), matrix(matrix), ins(matrix, 1) {}

	~scalar_assign()
	{
	    vampir_trace<3047> tracer;
	    if (with_comma) {
		MTL_CRASH_IF(r != num_rows(matrix), "Not all matrix entries initialized!");
	    } else {
		using std::min;
		if (src == math::zero(src)) // it is already set to zero
		    return;
		// Otherwise set diagonal (if square)
		for (SizeType i= 0, n= min(num_rows(matrix), num_cols(matrix)); i < n; i++)
		    ins[i][i] << src;
	    }
	}

	template <typename ValueSource>
	scalar_assign& operator, (ValueSource val)
	{
	    if (!with_comma) {
		with_comma= true;
		assert(r == 0 && c == 0);
		ins[r][c++] << src; // We haven't set v[0] yet
		if (c == num_cols(matrix)) 
		    c= 0, r++;
	    }
	    ins[r][c++] << val;
	    if (c == num_cols(matrix)) 
		c= 0, r++;
	    return *this;
	}

	Source         src;
	bool           with_comma;
	SizeType       r, c;
	Matrix&        matrix;
	inserter<Matrix> ins;
    };

    template <typename Source>
    typename boost::enable_if<boost::is_same<typename ashape::ashape<Source>::type, ashape::scal>, 
			      scalar_assign<Source> >::type
    operator=(Source src)
    {
	Matrix& matrix= static_cast<Matrix&>(*this);
	MTL_CRASH_IF(num_rows(matrix) * num_cols(matrix) == 0, 
		     "Trying to initialize a 0 by 0 matrix with a value!");
	set_to_zero(matrix);
	return scalar_assign<Source>(src, static_cast<Matrix&>(*this));
    }

    template <typename Source>
    Matrix& operator+=(const Source& src)
    {
	return density_plus_assign(src, mtl::traits::eval_dense< mat_mat_asgn_expr<Matrix, Source> >());
    }
    
    template <typename Source>
    Matrix& operator-=(const Source& src)
    {
	return density_minus_assign(src, mtl::traits::eval_dense< mat_mat_asgn_expr<Matrix, Source> >());
    }
    
    /// Scale matrix (in place) with scalar value or other matrix
    template <typename Factor>
    Matrix& operator*=(const Factor& alpha)
    {
	right_scale_inplace(static_cast<Matrix&>(*this), alpha);
	return static_cast<Matrix&>(*this);
    }

    // Element-wise scaling from right (i.e. like *= as elementwise)
    template <typename Factor>
    Matrix& ele_rscale(const Factor& alpha)
    {
	return density_ele_rscale(alpha, mtl::traits::eval_dense< mat_mat_asgn_expr<Matrix, Factor> >());
    }

    /// Divide matrix (in place) by scalar value
    // added by Hui Li
    template <typename Factor>
    Matrix& operator/=(const Factor& alpha)
    {
	divide_by_inplace(static_cast<Matrix&>(*this), alpha);
	return static_cast<Matrix&>(*this);
    }
};



template <typename Matrix, typename ValueType, typename SizeType>
struct const_crtp_matrix_bracket
{    
    template <typename T>
    typename boost::disable_if_c<boost::is_same<T, mtl::irange>::value || boost::is_same<T, mtl::iset>::value,
				 operations::bracket_proxy<Matrix, const Matrix&, ValueType> >::type
    operator[] (const T& row) const
    {
	return operations::bracket_proxy<Matrix, const Matrix&, ValueType>(static_cast<const Matrix&>(*this), row);
    }

    // Compiler error (later) if no sub_matrix function (or row vector resp.) available
    template <typename T>
    typename boost::enable_if<boost::is_same<T, mtl::irange>, operations::range_bracket_proxy<Matrix, const Matrix&, const Matrix> >::type
    operator[] (const T& row_range) const
    {
	return operations::range_bracket_proxy<Matrix, const Matrix&, const Matrix>(static_cast<const Matrix&>(*this), row_range);
    }

    operations::set_bracket_proxy<Matrix, const Matrix&, const Matrix>
    operator[] (const iset& row_set) const
    {
	return operations::set_bracket_proxy<Matrix, const Matrix&, const Matrix>(static_cast<const Matrix&>(*this), row_set);
    }
};

template <typename Matrix, typename ValueType, typename SizeType>
struct crtp_matrix_bracket 
{    
    operations::bracket_proxy<Matrix, const Matrix&, const ValueType&>
    operator[] (SizeType row) const
    {
        return operations::bracket_proxy<Matrix, const Matrix&, const ValueType&>(static_cast<const Matrix&>(*this), row);
    }

    template <typename T>
    typename boost::disable_if_c<boost::is_same<T, mtl::irange>::value || boost::is_same<T, mtl::iset>::value, 
			       operations::bracket_proxy<Matrix, Matrix&, ValueType&> >::type
    // operations::bracket_proxy<Matrix, Matrix&, ValueType&>
    operator[] (const T& row)
    {
        return operations::bracket_proxy<Matrix, Matrix&, ValueType&>(static_cast<Matrix&>(*this), row);
    }

    // Compiler error (later) if no sub_matrix function available
    operations::range_bracket_proxy<Matrix, const Matrix&, const Matrix>
    operator[] (const irange& row_range) const
    {
	return operations::range_bracket_proxy<Matrix, const Matrix&, const Matrix>(static_cast<const Matrix&>(*this), row_range);
    }

    // Compiler error (later) if no sub_matrix function available
    template <typename T>
    typename boost::enable_if<boost::is_same<T, mtl::irange>, operations::range_bracket_proxy<Matrix, Matrix&, Matrix> >::type
    // operations::range_bracket_proxy<Matrix, Matrix&, Matrix>
    operator[] (const T& row_range)
    {
	return operations::range_bracket_proxy<Matrix, Matrix&, Matrix>(static_cast<Matrix&>(*this), row_range);
    }

    operations::set_bracket_proxy<Matrix, const Matrix&, const Matrix>
    operator[] (const iset& row_set) const
    {
	return operations::set_bracket_proxy<Matrix, const Matrix&, const Matrix>(static_cast<const Matrix&>(*this), row_set);
    }
};

template <typename Matrix, typename ValueType, typename SizeType>
struct crtp_matrix_lvalue 
{ 
    // Function must be overwritten by Matrix if m(row, col) does not return a reference
    ValueType& lvalue(SizeType row, SizeType col)
    {
	return static_cast<Matrix&>(*this)(row, col);
    }   
};

template <typename Matrix, typename ValueType, typename SizeType>
struct const_crtp_base_matrix
  : public const_crtp_matrix_bracket<Matrix, ValueType, SizeType>
{};

template <typename Matrix, typename ValueType, typename SizeType>
struct mutable_crtp_base_matrix 
  : public crtp_matrix_bracket<Matrix, ValueType, SizeType>,
    public crtp_matrix_assign<Matrix, ValueType, SizeType>
{};

template <typename Matrix, typename ValueType, typename SizeType>
struct crtp_base_matrix 
  : boost::mpl::if_<boost::is_const<Matrix>,
		    const_crtp_base_matrix<Matrix, ValueType, SizeType>,
		    mutable_crtp_base_matrix<Matrix, ValueType, SizeType>
                   >::type
{};



}} // namespace mtl::matrix

#endif // MTL_CRTP_BASE_MATRIX_INCLUDE
