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

#ifndef MTL_TRIVIAL_INSERTER_INCLUDE
#define MTL_TRIVIAL_INSERTER_INCLUDE

#include <boost/numeric/mtl/operation/update.hpp>
#include <boost/numeric/mtl/operation/size.hpp>
#include <boost/numeric/mtl/matrix/element_matrix.hpp> 
#include <boost/numeric/mtl/matrix/element_array.hpp> 

namespace mtl { namespace detail {


// Matrix must have direct write access, i.e. matrix(row, col) must return a non-const reference
template <typename Matrix, typename Updater = mtl::operations::update_store<typename Matrix::value_type> >
struct trivial_inserter
{
    typedef trivial_inserter                            self;
    typedef Matrix                                      matrix_type;
    typedef typename matrix_type::size_type             size_type;
    typedef typename matrix_type::value_type            value_type;
    typedef operations::update_proxy<self, size_type>   proxy_type;
    
    explicit trivial_inserter(matrix_type& matrix, size_type) : matrix(matrix) {}

    proxy_type operator() (size_type row, size_type col)
    {
	return proxy_type(*this, row, col);
    }


  private:
    
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
	Updater() (matrix(row, col), val);
    }

    template <typename Modifier, typename Value>
    void modify(size_type row, size_type col, Value val)
    {
	Modifier() (matrix(row, col), val);
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

}} // namespace mtl::detail

#endif // MTL_TRIVIAL_INSERTER_INCLUDE
