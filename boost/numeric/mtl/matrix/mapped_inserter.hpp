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

#ifndef MTL_MATRIX_MAPPED_INSERTER_INCLUDE
#define MTL_MATRIX_MAPPED_INSERTER_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/update.hpp>

namespace mtl { namespace mat {

/// Inserter with shifted row and column indices
/** The main work is performed by the underlying base inserter whose type is given as template
    argument. **/
template <typename BaseInserter, typename Mapper > 
class mapped_inserter
{
  public:
    typedef mapped_inserter                                self;
    typedef typename BaseInserter::matrix_type              matrix_type;
    typedef typename Collection<matrix_type>::size_type     size_type;
    typedef operations::update_proxy<BaseInserter, size_type>   proxy_type;

    /// Constructor with matrix \p A, the mapping, and the slot size
    mapped_inserter(matrix_type& A, Mapper& map, size_type slot_size= 0)
      : ins(A, slot_size), map(map) {}

  private:
    struct bracket_proxy
    {
	bracket_proxy(BaseInserter& ref, Mapper& map, size_type row) 
          : ref(ref),map(map), row(row) {}
	
	proxy_type operator[](size_type col)
	{   return proxy_type(ref, row, map.col(col)); }

	BaseInserter&   ref;
        Mapper&         map;
	size_type row;
    };
    
  public:
    /// To be used in ins[r][c] << value;
    bracket_proxy operator[] (size_type row)   
    {	return bracket_proxy(ins, map, map.row(row));    }

    /// To be used in ins(r, c) << value;
    proxy_type operator() (size_type row, size_type col)  
    {	return proxy_type(ins, map.row(row),map.col(col));    }

    // update, modify and operator<< are used from BaseInserter

private:
    BaseInserter ins;
    Mapper&     map;
};

template< typename BaseInserter, typename Mapper, typename Elt, typename Parameters >
mapped_inserter< BaseInserter, Mapper >& operator<<(mapped_inserter< BaseInserter, Mapper >& minserter, 
    const compressed2D< Elt, Parameters >& rhs)
{

  using mtl::tag::major; using mtl::tag::nz; using mtl::begin; using mtl::end;
  namespace traits= mtl::traits;
  typedef compressed2D< Elt, Parameters >   Matrix;

  typename traits::row<Matrix>::type                                 row(rhs);
  typename traits::col<Matrix>::type                                 col(rhs);
  typename traits::const_value<Matrix>::type                         value(rhs);

  typedef typename traits::range_generator<major, Matrix>::type      cursor_type;
  typedef typename traits::range_generator<nz, cursor_type>::type    icursor_type;

  for (cursor_type cursor = begin<major>(rhs), cend = end<major>(rhs); cursor != cend; ++cursor)
    for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor); icursor != icend; ++icursor)
      minserter[row(*icursor)][col(*icursor)]<< value(*icursor);

  return minserter;
}

} // namespace matrix

} // namespace mtl

#endif // MTL_MATRIX_MAPPED_INSERTER_INCLUDE
