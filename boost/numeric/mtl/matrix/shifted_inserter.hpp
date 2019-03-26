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

#ifndef MTL_MATRIX_SHIFTED_INSERTER_INCLUDE
#define MTL_MATRIX_SHIFTED_INSERTER_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/update.hpp>

namespace mtl { namespace mat {

/// Inserter with shifted row and column indices
/** The main work is performed by the underlying base inserter whose type is given as template
    argument. **/
template <typename BaseInserter> 
class shifted_inserter
{
  public:
    typedef shifted_inserter                                self;
    typedef typename BaseInserter::matrix_type              matrix_type;
    typedef typename Collection<matrix_type>::size_type     size_type;
    typedef operations::update_proxy<BaseInserter, size_type>   proxy_type;

    /// Constructor with matrix \p A, the slot size and the offsets
    shifted_inserter(matrix_type& A, size_type slot_size= 0,
		     size_type row_offset= 0, size_type col_offset= 0)
      : ins(A, slot_size), row_offset(row_offset), col_offset(col_offset) {}

    void set_row_offset(size_type ro) { row_offset= ro; } ///< Change row offset
    void set_col_offset(size_type co) { col_offset= co; } ///< Change column offset

    size_type get_row_offset() const { return row_offset; } ///< Get row offset
    size_type get_col_offset() const { return col_offset; } ///< Get column offset

  private:
    struct bracket_proxy
    {
	bracket_proxy(BaseInserter& ref, size_type row, size_type col_offset) : ref(ref), row(row), col_offset(col_offset) {}
	
	proxy_type operator[](size_type col)
	{   return proxy_type(ref, row, col+col_offset); }

	BaseInserter&   ref;
	size_type       row, col_offset;
    };
    
  public:
    /// To be used in ins[r][c] << value;
    bracket_proxy operator[] (size_type row)      
    {	return bracket_proxy(ins, row+row_offset, col_offset);    } 

    /// To be used in ins(r, c) << value;
    proxy_type operator() (size_type row, size_type col)  
    {	return proxy_type(ins, row+row_offset, col+col_offset);    }

    // update, modify and operator<< are used from BaseInserter

  private:
    BaseInserter ins;
    size_type    row_offset, col_offset;
};

}} // namespace mtl::matrix

#endif // MTL_MATRIX_SHIFTED_INSERTER_INCLUDE
