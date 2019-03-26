
#ifndef MTL_DIST_BLOCK_DIAGONAL2D_INCLUDE
#define MTL_DIST_BLOCK_DIAGONAL2D_INCLUDE



#include <vector>
#include <cassert>
#include <limits>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/bool.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/assert.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/mtl/matrix/mat_expr.hpp>


namespace mtl {

  namespace mat {
    

/// Block diagonal matrix structure
/** Blocks can be any existing matrix_type in mtl4. **/
template <typename Matrix>
class block_diagonal2D 
  : public mat_expr< block_diagonal2D<Matrix> >
{
  public:
 
    typedef unsigned int	      size_type;
    typedef Matrix		      block_type;
    typedef std::vector<Matrix>       block_matrix_type;
    typedef std::vector<size_type>    index_vector_type;
    typedef block_diagonal2D<Matrix>  self;
    typedef typename Collection<Matrix>::value_type value_type;

    /// Constructor: number or rows and columns and optionally estimated number of blocks
    explicit block_diagonal2D(size_type rows, size_type cols, size_type init_size= 0)
      : nrows(rows), ncols(cols), nb_blocks(0), my_nnz(0),
	min_ind(std::numeric_limits<size_type>::max()), max_ind(std::numeric_limits<size_type>::min()),
	compact_heap(0)
    {
	start_block.reserve(init_size); 
	end_block.reserve(init_size); 
	blocks.reserve(init_size);
    }

    ~block_diagonal2D() { delete[] compact_heap; }

    /// Returns the global number of rows
    size_type num_rows() const {
	return nrows ;
    }
    /// Returns the global number of coluums
    size_type num_cols() const {
	return ncols ;
    }
    /// Returns the global number of blocks
    size_type num_blocks() const {
	return  nb_blocks; 
    }

    /// Minimal index
    size_type min_index() const { return min_ind; }

    /// Maximal index
    size_type max_index() const { return max_ind; }


    block_type const& block(size_type i) const
    {
	MTL_CRASH_IF(is_negative(i) || i >= nb_blocks, "Index out of range!");
	return blocks[i];
    }

    /// Insert a block from start x start to end x end
    void insert(size_type start, size_type end, const block_type& A) 
    {
	MTL_CRASH_IF(start > end, "Logic error");
	MTL_CRASH_IF(is_negative(start) || end > nrows || end > ncols, "Index out of range!");
	assert(compact_heap == 0); // insertion after make_compact; might be relaxed later
    
	start_block.push_back(start);
	end_block.push_back(end);
	blocks.push_back(A); 
	my_nnz+= size_type(A.nnz());
	nb_blocks++;
	
	if (start < min_ind)
	    min_ind= start;
	if (end > max_ind)
	    max_ind= end;
#if 0
	boost::mpi::communicator world;
	if (world.rank() == 0) {
	    std::cout << "block_diag.insert " << nb_blocks-1 << "th block: start == "
		      << start_block << ", end == " << end_block << "block is:\n"
		      << blocks.back();
	}
#endif
    } 

    /// Element A[i][j] by summing over all blocks, use only for debugging because it is very slow
    value_type operator()(size_type i, size_type j) const
    {
	value_type s= value_type(0);
	for (std::size_t b= 0; b < blocks.size(); ++b) {
	    size_type st= start_block[b], e= end_block[b];
	    if (st <= i && st <= j && i < e && j < e)
		s+= blocks[b][i - st][j - st];
	}
	return s;
    }
 
    /// Memory of inserted 
    void make_compact(boost::mpl::true_)
    {
	assert(compact_heap == 0); // might be relaxed later
	
	size_type entries= 0;
	for (size_type i= 0; i < nb_blocks; i++)
	    entries+= size(blocks[i]);
	compact_heap= new value_type[entries];

	size_type pos= 0;
	for (size_type i= 0; i < nb_blocks; i++) {
	    size_type s= end_block[i] - start_block[i];
	    block_type tmp(s, s, compact_heap + pos);
	    tmp= blocks[i];
	    pos+= size(blocks[i]);
	    swap(blocks[i], tmp);
	}
	assert(pos == entries);
    }
    void make_compact(boost::mpl::false_) {}

    void make_compact() { make_compact(boost::is_same<typename mtl::traits::category<block_type>::type, tag::dense2D>()); }
 
    /// Number of non-zeros (accumulated over blocks)
    size_type nnz() const { return my_nnz; }

    template <typename VectorIn, typename VectorOut>
    void add_mult(const VectorIn& x, VectorOut& y) const 
    {
	mtl::vampir_trace<3062> tracer;
	MTL_CRASH_IF(ncols != size(x) || nrows != size(y), "Incompatible size!"); 

	// set_to_zero(y);
	for(size_type i= 0; i < nb_blocks; i++) {
	    mtl::irange r(start_block[i], end_block[i]);
	    y[r]+= blocks[i] * x[r];
	}
    }

    ///block_diagonal-matrix  times cvec
    template <typename VectorIn, typename VectorOut>
    void mult(const VectorIn& x, VectorOut& y) const 
    {
	mtl::vampir_trace<3062> tracer;
	MTL_CRASH_IF(ncols != size(x) || nrows != size(y), "Incompatible size!"); 

	set_to_zero(y);
	for(size_type i= 0; i < nb_blocks; i++) {
	    mtl::irange r(start_block[i], end_block[i]);
	    y[r]= blocks[i] * x[r];
	}
    }

    template <typename VectorIn>
    struct multiplier
      : mtl::vec::assigner<multiplier<VectorIn> >
    {
	explicit multiplier(const self& P, const VectorIn& x) : P(P), x(x) {}

	template <typename VectorOut>
	void assign_to(VectorOut& y) const
	{   P.mult(x, y); }
	
	const self& P;
	const VectorIn& x;
    };
  
 
    template <typename VectorIn>
    multiplier<VectorIn> operator*(const VectorIn& x) // const
    {  return multiplier<VectorIn>(*this, x); }


  private:
    size_type           nrows, ncols, nb_blocks, my_nnz, min_ind, max_ind;
    index_vector_type  	start_block, end_block;
    block_matrix_type   blocks;
    value_type*         compact_heap; 
};

// ================
// Free functions
// ================


/// Number of rows
template <typename Value>
typename block_diagonal2D<Value>::size_type
inline num_rows(const block_diagonal2D<Value>& matrix)
{
    return matrix.num_rows();
}

/// Number of columns
template <typename Value>
typename block_diagonal2D<Value>::size_type
inline num_cols(const block_diagonal2D<Value>& matrix)
{
    return matrix.num_cols();
}

/// Size of the matrix, i.e. the number of row times columns
template <typename Value>
typename block_diagonal2D<Value>::size_type
inline size(const block_diagonal2D<Value>& matrix)
{
    return matrix.num_cols() * matrix.num_rows();
}

// /// Distributed Vector of start block entrys in the matrix
// template <typename Value>
// typename block_diagonal2D<Value>::index_type
// inline start(const block_diagonal2D<Value>& matrix)
// {
//     return matrix.start();
// }

// /// Distributed Vector of en block entrys in the matrix
// template <typename Value>
// typename block_diagonal2D<Value>::index_type
// inline end(const block_diagonal2D<Value>& matrix)
// {
//     return matrix.end();
// }




} // namespace matrix
} // namespace mtl

#endif // MTL_DIST_BLOCK_DIAGONAL2D_INCLUDE

