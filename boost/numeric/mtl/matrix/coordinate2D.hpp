
#ifndef MTL_COORDINATE2D_INCLUDE
#define MTL_COORDINATE2D_INCLUDE

#include <vector>
#include <cassert>
#include <boost/mpl/bool.hpp>


#include <boost/numeric/mtl/operation/is_negative.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/operation/sort.hpp>
#include <boost/numeric/mtl/operation/iota.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>

namespace mtl {  namespace mat {
    
  
/// Sparse matrix structure in coordinate format
template <typename Value, typename Parameters = mat::parameters<> >
class coordinate2D 
  : public base_matrix<Value, Parameters>,
    public const_crtp_base_matrix< coordinate2D<Value, Parameters>, Value, typename Parameters::size_type >,
    public crtp_matrix_assign< coordinate2D<Value, Parameters>, Value, typename Parameters::size_type >,
    public mat_expr< coordinate2D<Value, Parameters> >
{
  public:
    typedef Value                                                       value_type;
    typedef Value&                                                      reference;
    typedef Value const&                                                const_reference;
    typedef typename Parameters::size_type				  size_type;
    typedef typename Parameters::dimensions                             dim_type;
    typedef typename Parameters::orientation                            orientation;      

    typedef coordinate2D                                                self;
    typedef base_matrix<Value, Parameters>                              super;
    typedef crtp_matrix_assign<self, Value, size_type>                  assign_base;

    typedef std::vector< size_type >                                    row_index_array_type ;
    typedef std::vector< size_type >                                    column_index_array_type ;
    typedef std::vector< value_type >                                   value_array_type ;

    /// Common constructor
    explicit coordinate2D(size_type nrows, size_type ncols, size_type expected= 0)
      : super(dim_type(nrows, ncols))
    {
	if (expected > 0) {
	    rows.reserve(expected);
	    cols.reserve(expected);
	    values.reserve(expected);
	}
	my_is_sorted= true; 
    } 
    
    using assign_base::operator=;    
  
    size_type nnz() const { return rows.size(); } ///< Number of non-zeros

    value_array_type const& value_array() const { return values; } ///< Array of values (const)
    value_array_type& value_array() { return values; }             ///< Array of values (mutable)

    row_index_array_type const& row_index_array() const { return rows; }       ///< Array of rows  (const)
    column_index_array_type const& column_index_array() const {	return cols; } ///< Array of columns (const)
  
    row_index_array_type& row_index_array() { return rows; }       ///< Array of rows   (mutable)
    column_index_array_type& column_index_array() { return cols; } ///< Array of columns  (mutable)

    /// Drop all entries
    void make_empty()
    {
        rows.resize(0); cols.resize(0); values.resize(0);
        my_is_sorted= true; // haha
    }
    
    /// Insert an entry at the end of the row-,col- and value-array 
    void push_back(size_type r, size_type c, const_reference v) 
    {
	rows.push_back(r); cols.push_back(c); values.push_back(v); my_is_sorted= false;
    } 
  
    /// Insert an entry at the end of the row-,col- and value-array, like push_back
    void insert(size_type r, size_type c, const_reference v) {	push_back(r, c, v); }
  
    /// Whether the entries are sorted
    bool is_sorted() const { return my_is_sorted; }

    /// sorting standard by rows
    void sort() 
    { 
	if (nnz() > 0) 
	    sort(mtl::traits::is_row_major<Parameters>());  
	my_is_sorted= true;
    }

  private:
    // sorting by rows
    void sort(boost::mpl::true_)
    {  
	mtl::vec::sort_xy(rows, cols, values);
    }

    // sorting by columns
    void sort(boost::mpl::false_)
    {
	mtl::vec::sort_xy(cols, rows, values);
    }
  
    template <typename OStream, typename Vector>
    void print_stl_vector(OStream& os, const Vector& v) const
    {
	os << "[";
	for (unsigned i= 0; i < v.size(); i++)
	    os << v[i] << (i+1 < v.size() ? "," : "");
	os << "]\n";
    }

  public:
    template <typename OStream>
    void print_internal(OStream& os) const
    {
       	os << "rows   = "; print_stl_vector(os, rows);
	os << "cols   = "; print_stl_vector(os, cols);
	os << "values = "; print_stl_vector(os, values);
    }

    void print_internal() const { print(std::cout); }

    ///operator * for  vector= coordinaten-matrix * vector
    template <typename Vector >
    Vector operator*(const Vector& x)
    {
	
	Vector res(this->num_rows());
	res= 0;
	for (size_type i= 0; i < nnz(); i++)
	    res[rows[i]]+=  values[i] * x[cols[i]];
	return res;
    }
  
    value_type operator() (const size_type r, const size_type c) const
    {
	MTL_CRASH_IF(is_negative(r) || r >= this->num_rows() 
		  || is_negative(c) || c >= this->num_cols(), "Index out of range!");

#if 0
	if (my_is_sorted)
	    return find(r, c, mtl::traits::is_row_major<Parameters>());
#endif

	for (size_type i= 0; i < nnz(); i++) 
	    if (rows[i] == r && cols[i] == c)
		return values[i];
	return value_type(0);
    }

    template <typename Updater>
    void compress(Updater up)
    {
	if (!my_is_sorted)
	    sort();

	size_type i= 0, j= 1, end= rows.size();
	for (; j < end; ++j) 
	    if (rows[i] == rows[j] && cols[i] == cols[j]) {
		up(values[i], values[j]);
	    } else {
		i++;
		if (i != j) {
		    rows[i]= rows[j];
		    cols[i]= cols[j];
		    values[i]= values[j];
		}
	    }
	if (end > 0) i++;
	rows.resize(i);
	cols.resize(i);
	values.resize(i);
   }

    template <typename Matrix, typename Updater> friend struct coordinate2D_inserter;

  private:
    row_index_array_type      rows;
    column_index_array_type   cols;
    value_array_type          values;
    bool                      my_is_sorted;
};

// ================
// Free functions
// ================


/// Number of rows
template <typename Value, typename Parameters>
typename coordinate2D<Value, Parameters>::size_type
inline num_rows(const coordinate2D<Value, Parameters>& matrix)
{
    return matrix.num_rows();
}

/// Number of columns
template <typename Value, typename Parameters>
typename coordinate2D<Value, Parameters>::size_type
inline num_cols(const coordinate2D<Value, Parameters>& matrix)
{
    return matrix.num_cols();
}

/// Size of the matrix, i.e. the number of row times columns
template <typename Value, typename Parameters>
typename coordinate2D<Value, Parameters>::size_type
inline size(const coordinate2D<Value, Parameters>& matrix)
{
    return matrix.num_cols() * matrix.num_rows();
}

/// Number of NoZeros of the matrix
template <typename Value, typename Parameters>
typename coordinate2D<Value, Parameters>::size_type
inline nnz(const coordinate2D<Value, Parameters>& matrix)
{
    return matrix.nnz();
}


template <typename Matrix, 
	  typename Updater = mtl::operations::update_store<typename Matrix::value_type> >
struct coordinate2D_inserter
{
    typedef coordinate2D_inserter                       self;
    typedef Matrix                                      matrix_type;
    typedef typename matrix_type::size_type             size_type;
    typedef typename matrix_type::value_type            value_type;
    typedef operations::update_proxy<self, size_type>   proxy_type;
    
    // We only support storing so far !!!
    // STATIC_ASSERT((boost::is_same<Updater, mtl::operations::update_store<value_type> >::value), "We only support storing so far");

    coordinate2D_inserter(matrix_type& matrix, size_type slot_size= 1) 
      : matrix(matrix) 
    {
	std::size_t ns= slot_size * matrix.dim1();
	if (ns > matrix.nnz()) {
	    matrix.rows.reserve(ns);
	    matrix.cols.reserve(ns);
	    matrix.values.reserve(ns);
	}
    }

    ~coordinate2D_inserter() { matrix.compress(Updater()); }
    
 private:

    struct update_proxy
    {
	// self is type of inserter not update_proxy !!!
	update_proxy(self& ref, size_type row, size_type col) : ref(ref), row(row), col(col) {}

	template <typename Value>
	update_proxy& operator<< (Value const& val)
	{
	    ref.matrix.push_back(row, col, val);
	    return *this;
	}
	self& ref;
	size_type row, col;
    };
    
    proxy_type operator() (size_type row, size_type col)
    {
	return proxy_type(*this, row, col);
    }

    
    struct bracket_proxy
    {
	bracket_proxy(self& ref, size_type row) : ref(ref), row(row) {}
	
	proxy_type operator[](size_type col)
	{
	    return proxy_type(ref, row, col);
	}

	self&      ref;
	size_type  row;
    };

  public:

    bracket_proxy operator[] (size_type row)
    {
	return bracket_proxy(*this, row);
    }

    template <typename Value>
    void update(size_type row, size_type col, Value val)
    {
	matrix.push_back(row, col, val);
    }

    template <typename Modifier, typename Value>
    void modify(size_type row, size_type col, Value val)
    {
	matrix.push_back(row, col, val);
    }

    template <typename EMatrix, typename Rows, typename Cols>
    self& operator<< (const mat::element_matrix_t<EMatrix, Rows, Cols>& elements)
    {
	using mtl::size;
	for (unsigned ri= 0; ri < size(elements.rows); ri++)
	    for (unsigned ci= 0; ci < size(elements.cols); ci++)
		update (elements.rows[ri], elements.cols[ci], elements.matrix(ri, ci));
	return *this;
    }

    template <typename EMatrix, typename Rows, typename Cols>
    self& operator<< (const mat::element_array_t<EMatrix, Rows, Cols>& elements)
    {
	using mtl::size;
	for (unsigned ri= 0; ri < size(elements.rows); ri++)
	    for (unsigned ci= 0; ci < size(elements.cols); ci++)
		update (elements.rows[ri], elements.cols[ci], elements.array[ri][ci]);
	return *this;
    }

  protected:
    matrix_type&         matrix;
};

struct coordinate_key
{
    typedef std::size_t                               size_t;

    explicit coordinate_key(size_t offset) : offset(offset) {}

    bool operator== (coordinate_key const& other) const { return offset == other.offset; }
    bool operator!= (coordinate_key const& other) const { return offset != other.offset; }
    
    size_t offset;    
};


// Cursor over every element
template <typename Value, typename Parameters>
struct coordinate_minor_cursor 
 : public coordinate_key 
{
    typedef coordinate_minor_cursor<Value, Parameters>   self;
    typedef typename Parameters::size_type               size_type;
    typedef const coordinate2D<Value, Parameters>&       matrix_ref_type;
    static const int                                     level= 2;

    coordinate_minor_cursor(matrix_ref_type ref, size_type offset) 
      : coordinate_key(offset), ref(ref)  {}

    bool operator!=(const self& that) const
    {
	assert(&ref == &that.ref);
	return this->offset != that.offset;
    }

    self& operator++() { this->offset++; return *this; }
    coordinate_key operator*() const { return *this; }

    matrix_ref_type ref;
};


template <typename Value, typename Parameters>
struct coordinate_major_cursor 
{
    typedef coordinate_major_cursor<Value, Parameters>   self;
    typedef typename Parameters::size_type               size_type;
    typedef const coordinate2D<Value, Parameters>&       matrix_ref_type;
    typedef coordinate_minor_cursor<Value, Parameters>   inner_cursor;
    static const int                                     level= 2;

    void find_next_offset(boost::mpl::true_)
    {
	size_type i= offset;
	for ( ; i < nnz(ref) && ref.row_index_array()[i] <= major; i++) ;
	next_offset= i;
    }

    void find_next_offset(boost::mpl::false_)
    {
	size_type i= offset;
	for ( ; i < nnz(ref) && ref.col_index_array()[i] <= major; i++) ;
	next_offset= i;
    }

    void find_next_offset() { find_next_offset(mtl::traits::is_row_major<Parameters>()); }

    coordinate_major_cursor(matrix_ref_type ref, size_type major, size_type offset) 
      : ref(ref), major(major), offset(offset)
    {
	find_next_offset();
    }

    bool operator!=(const self& that) const
    {
	assert(&ref == &that.ref);
	return this->offset != that.offset;
    }

    self& operator++() 
    { 
	offset= next_offset; 
	major++;
	find_next_offset();
	return *this; 
    }

    matrix_ref_type ref;
    size_type       major, offset, next_offset;
};

template <typename Value, typename Parameters>
struct coordinate_minor_range_generator
{
    typedef coordinate_major_cursor<Value, Parameters>   outer_cursor_type;
    typedef coordinate_minor_cursor<Value, Parameters>   type;
    static const int                                     level= 2;

    type begin(outer_cursor_type c) const { return type(c.ref, c.offset); }
    type end(outer_cursor_type c) const { return type(c.ref, c.next_offset); }
};

template <typename Value, typename Parameters>
struct coordinate_row_range_generator
{
    typedef const coordinate2D<Value, Parameters>&       matrix_ref_type;
    typedef coordinate_major_cursor<Value, Parameters>   type;
    typedef complexity_classes::linear_cached            complexity;
    static const int                                     level= 1;

    type begin(matrix_ref_type A) const { return type(A, 0, 0); }
    type end(matrix_ref_type A) const { return type(A, num_rows(A), nnz(A)); }
};

template <typename Value, typename Parameters>
struct coordinate_col_range_generator
{
    typedef const coordinate2D<Value, Parameters>&       matrix_ref_type;
    typedef coordinate_major_cursor<Value, Parameters>   type;
    typedef complexity_classes::linear_cached            complexity;
    static const int                                     level= 1;

    type begin(matrix_ref_type A) const { return type(A, 0, 0); }
    type end(matrix_ref_type A) const { return type(A, num_cols(A), nnz(A)); }
};


}} // namespace mtl::matrix


namespace mtl { namespace traits {
	
    // Cursor over all rows
    // Supported if row major matrix
    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::row, mat::coordinate2D<Value, Parameters> >
      : boost::mpl::if_<
	    boost::is_same<typename Parameters::orientation, row_major>
	  , mat::coordinate_row_range_generator<Value, Parameters>
 	  , range_generator<tag::unsupported, mat::coordinate2D<Value, Parameters> >
        >::type {};	

    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::col, mat::coordinate2D<Value, Parameters> >
      : boost::mpl::if_<
	    boost::is_same<typename Parameters::orientation, col_major>
	  , mat::coordinate_col_range_generator<Value, Parameters>
 	  , range_generator<tag::unsupported, mat::coordinate2D<Value, Parameters> >
        >::type {};	

    template <class Value, class Parameters>
    struct range_generator<glas::tag::nz, mat::coordinate_major_cursor<Value, Parameters> >
      : mat::coordinate_minor_range_generator<Value, Parameters>
    {};


}} // namespace mtl::traits

namespace mtl {
	using mat::coordinate2D;
}

#endif // MTL_COORDINATE2D_INCLUDE

