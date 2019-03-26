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

#ifndef MTL_MATRIX_SPARSE_BANDED_INCLUDE
#define MTL_MATRIX_SPARSE_BANDED_INCLUDE

#include <algorithm>
#include <cassert>
#include <ostream>
#include <vector>
#include <boost/type_traits/make_signed.hpp>

#include <boost/numeric/mtl/matrix/dimension.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/base_matrix.hpp>
#include <boost/numeric/mtl/matrix/crtp_base_matrix.hpp>
#include <boost/numeric/mtl/matrix/mat_expr.hpp>
#include <boost/numeric/mtl/operation/is_negative.hpp>
#include <boost/numeric/mtl/operation/update.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/maybe.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>

namespace mtl { namespace mat {

/// Sparse banded matrix class
template <typename Value, typename Parameters = mat::parameters<> >
class sparse_banded
  : public base_matrix<Value, Parameters>,
    public const_crtp_base_matrix< sparse_banded<Value, Parameters>, Value, typename Parameters::size_type >,
    public crtp_matrix_assign< sparse_banded<Value, Parameters>, Value, typename Parameters::size_type >,
    public mat_expr< sparse_banded<Value, Parameters> >
{
    MTL_STATIC_ASSERT((mtl::traits::is_row_major<Parameters>::value), "Only row-major sparse banded matrices supported so far.");

    typedef std::size_t                                size_t;
    typedef base_matrix<Value, Parameters>             super;
    typedef sparse_banded<Value, Parameters>           self;

  public:
    typedef Value                                      value_type;
    typedef typename Parameters::size_type             size_type;
    typedef typename boost::make_signed<size_type>::type  band_size_type;
    typedef crtp_matrix_assign<self, Value, size_type> assign_base;

    /// Construct matrix of dimension \p nr by \p nc
    explicit sparse_banded(size_type nr= 0, size_type nc= 0) 
      : super(non_fixed::dimensions(nr, nc)), data(0), inserting(false)
    {}

    /// Copy from other types
    template <typename MatrixSrc>
    explicit sparse_banded(const MatrixSrc& src) : data(0), inserting(false)
    {	*this= src;    }


    ~sparse_banded() { delete[] data; }
    void check() const { MTL_DEBUG_THROW_IF(inserting, access_during_insertion()); }
    void check(size_type MTL_DEBUG_ARG(r), size_type MTL_DEBUG_ARG(c)) const
    {
	check();
	MTL_DEBUG_THROW_IF(is_negative(r) || r >= this->num_rows() 
			   || is_negative(c) || c >= this->num_cols(), index_out_of_range());
    }

    using assign_base::operator=;

    void make_empty() ///< Delete all entries
    {
	delete[] data; data= 0;
	bands.resize(0);
    }

    /// Change dimension to \p r by \p c; deletes all entries
    void change_dim(size_type r, size_type c) 
    {
	make_empty();
	super::change_dim(r, c);
    }

    /// Offset of entry r,c if in matrix (otherwise offset of next entry)
    utilities::maybe<size_type> offset(size_type r, size_type c) const
    {
	check(r, c);
	band_size_type dia= band_size_type(c) - band_size_type(r);
	size_type      b= find_dia(dia);
	return utilities::maybe<size_type>(r * bands.size() + b, size_t(b) < bands.size() && bands[b] == dia);
    }

    /// Value of matrix entry
    value_type operator()(size_type r, size_type c) const
    {
	using math::zero;
	utilities::maybe<size_type> o= offset(r, c);
	return o ? data[o.value()] : zero(value_type());
    }

    /// L-value reference of stored matrix entry
    /** To be used with care; in debug mode, exception is thrown if entry is not found **/
    value_type& lvalue(size_type r, size_type c) 
    {
	utilities::maybe<size_type> o= offset(r, c);
	MTL_DEBUG_THROW_IF(!o, runtime_error("Entry doesn't exist."));
	return data[o.value()];
    }

    /// Print matrix on \p os
    friend std::ostream& operator<<(std::ostream& os, const self& A)
    {
	const size_type bs= A.bands.size(), nc= A.num_cols(), s= bs * A.num_rows();
	print_size_max p;
	for (size_type i= 0; i < s; i++)
	    p(A.data[i]);

	for (size_type r= 0; r < A.num_rows(); r++) {
	    os << '[';
	    size_type b= 0;
	    while (b < bs && -A.bands[b] > long(r)) b++;
	    for (size_type c= 0; c < nc; c++) {
		os.width(p.max);
		if (b == bs || long(c) - long(r) < A.bands[b])
		    os << Value(0);
		else
		    os << A.data[r * bs + b++];
		os << (c + 1 == nc ? "]\n" : " ");
	    }
	}
	return os;
    }

    /// Number of structural non-zeros (i.e. stored entries) which is the the number of bands times rows
    size_type nnz() const { return bands.size() * this->num_rows(); }
    size_type find_major(size_type offset) const { return offset % bands.size(); }

    friend size_t num_rows(const self& A) { return A.num_rows(); }
    friend size_t num_cols(const self& A) { return A.num_cols(); }

    std::vector<band_size_type> const& ref_bands() const { return bands; } ///< Refer bands vector [advanced]
    const value_type*                  ref_data() const { return data; } ///< Refer data pointer [advanced]

  private:
    size_type find_dia(band_size_type dia) const
    {
	size_type i= 0;
	while (i < size_type(bands.size()) && dia > bands[i]) i++;
	return i;
    }

    template <typename, typename, typename> friend struct sparse_banded_inserter;

    std::vector<band_size_type>    bands;
    value_type*                    data;
    bool                           inserting;
};

/// Inserter for sparse banded matrix
template <typename Value, typename Parameters, typename Updater = mtl::operations::update_store<Value> >
struct sparse_banded_inserter
{
    typedef Value                                                            value_type;
    typedef typename Parameters::size_type                                   size_type;
    typedef sparse_banded<Value, Parameters>                                 matrix_type;
    typedef sparse_banded_inserter                                           self;
    typedef typename matrix_type::band_size_type                             band_size_type;
    typedef operations::update_proxy<self, size_type>                        proxy_type;

  private:
    struct bracket_proxy
    {
	bracket_proxy(self& ref, size_type row) : ref(ref), row(row) {}
	
	proxy_type operator[](size_type col) { return proxy_type(ref, row, col); }

	self&      ref;
	size_type  row;
    };

    void insert_dia(size_type dia_band, band_size_type dia)
    {
	using std::swap;
	// empty vector and entry at the end
	diagonals.push_back(std::vector<Value>(num_rows(A), Value(0)));
	A.bands.push_back(dia);
	
	for (size_type i= diagonals.size() - 1; i > dia_band; i--) {
	    swap(diagonals[i-1], diagonals[i]);
	    swap(A.bands[i-1], A.bands[i]);
	}
    }

  public:
    /// Construct inserter for matrix \p A; second argument for slot_size ignored
    sparse_banded_inserter(matrix_type& A, size_type) : A(A) 
    {
	MTL_THROW_IF(A.inserting, runtime_error("Two inserters on same matrix"));
	MTL_THROW_IF(A.data, runtime_error("Inserter does not support modifications yet."));
	A.inserting= true;
    }

    ~sparse_banded_inserter()
    {
	const size_type bs= A.bands.size();
	Value* p= A.data= new Value[A.bands.size() * num_rows(A)];
	for (size_type r= 0; r < num_rows(A); r++)
	    for (size_type b= 0; b < bs; b++)
		*p++= diagonals[b][r];
	assert(p - A.data == long(A.bands.size() * num_rows(A)));
	A.inserting= false;
    }

    /// Proxy to insert into A[row][col]
    bracket_proxy operator[] (size_type row)
    {
	return bracket_proxy(*this, row);
    }

    /// Proxy to insert into A[row][col]
    proxy_type operator() (size_type row, size_type col)
    {
	return proxy_type(*this, row, col);
    }

    /// Modify A[row][col] with \p val using \p Modifier
    template <typename Modifier>
    void modify(size_type row, size_type col, Value val)
    {
	MTL_DEBUG_THROW_IF(is_negative(row) || row >= num_rows(A) || is_negative(col) || col >= num_cols(A), 
			   index_out_of_range());
	const band_size_type dia= col - row;
	const size_type      dia_band= A.find_dia(dia);

	if (dia_band == size_type(A.bands.size()) || dia != A.bands[dia_band])
	    insert_dia(dia_band, dia);

	Modifier()(diagonals[dia_band][row], val);
    }

    /// Modify A[row][col] with \p val using the class' updater
    void update(size_type row, size_type col, Value val)
    {
	using math::zero;
	modify<Updater>(row, col, val);
    }

  private:
    matrix_type&                      A;
    std::vector<std::vector<Value> >  diagonals;
};


template <typename Value, typename Parameters>
struct sparse_banded_key
{
    typedef std::size_t                                   size_t;
    typedef sparse_banded_key                             self;
    typedef mtl::mat::sparse_banded<Value, Parameters> matrix_type;
    
    explicit sparse_banded_key(const matrix_type& A, size_t offset) 
      : A(A), offset(offset) {}

    bool operator== (sparse_banded_key const& other) const { return offset == other.offset; }
    bool operator!= (sparse_banded_key const& other) const { return !(*this == other); }

    const matrix_type&   A;
    size_t               offset;
};

// Cursor over every element
template <typename Value, typename Parameters>
struct sparse_banded_cursor 
  : public sparse_banded_key<Value, Parameters>
{
    typedef sparse_banded_key<Value, Parameters>          base;
    typedef mtl::mat::sparse_banded<Value, Parameters> matrix_type;
    typedef sparse_banded_cursor                          self;

    explicit sparse_banded_cursor(const matrix_type& A, size_t offset) 
      : base(A, offset)    {}

    self& operator++() { ++this->offset; return *this;  }
    const base& operator* () const { return *this; }
};

// ================
// Free functions
// ================

template <typename Value, typename Parameters>
inline std::size_t size(const sparse_banded<Value, Parameters>& matrix)
{
    return std::size_t(matrix.num_cols()) * std::size_t(matrix.num_rows());
}

}} // namespace mtl::matrix

namespace mtl {
    using mat::sparse_banded;
}


// ================
// Range generators
// ================

namespace mtl { namespace traits {

    // Cursor over all rows
    // Supported if row major matrix
    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::row, sparse_banded<Value, Parameters> >
      : detail::all_rows_range_generator<sparse_banded<Value, Parameters>, complexity_classes::linear_cached>
    {};

    template <class Value, class Parameters>
    struct range_generator<glas::tag::nz, 
			   detail::sub_matrix_cursor<sparse_banded<Value, Parameters>, glas::tag::row, 2> >
    {
	typedef typename Parameters::size_type                                                 size_type;
	typedef detail::sub_matrix_cursor<sparse_banded<Value, Parameters>, glas::tag::row, 2> cursor_type;
	typedef complexity_classes::linear_cached                                              complexity;
	typedef mtl::mat::sparse_banded_cursor<Value, Parameters>                           type;
	static int const                                                                       level = 1;
	MTL_STATIC_ASSERT((is_row_major<Parameters>::value), "Only row-major sparse banded matrices supported so far.");

	type begin(cursor_type const& cursor) const
	{
	    size_type r= cursor.key, b= 0, bs= cursor.ref.ref_bands().size();
	    while (b < bs && -cursor.ref.ref_bands()[b] > long(r)) b++;
	    return type(cursor.ref, r * bs + b);
	}
	type end(cursor_type const& cursor) const
	{
	    size_type r= cursor.key, bs= cursor.ref.ref_bands().size(), b= bs, nc= cursor.ref.num_cols();
	    while (b > 0 && cursor.ref.ref_bands()[b-1] + long(r) >= long(nc)) b--;
	    return type(cursor.ref, r * bs + b);
	}

	// type lower_bound(cursor_type const& cursor, size_type position) const {}
    };


    // Cursor over all rows or columns, depending which one is major
    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::major, sparse_banded<Value, Parameters> >
      : range_generator<glas::tag::row, sparse_banded<Value, Parameters> >
    {};

}} // namespace mtl::traits

#endif // MTL_MATRIX_SPARSE_BANDED_INCLUDE
