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

#ifndef MTL_VECTOR_MAPPED_INSERTER_INCLUDE
#define MTL_VECTOR_MAPPED_INSERTER_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/update.hpp>

namespace mtl {
namespace vec {
  
/// Inserter with shifted indices
/** The main work is performed by the underlying base inserter whose type is given as template
    argument. **/
template <typename BaseInserter, typename Mapper > 
class mapped_inserter
{
  public:
    typedef mapped_inserter                                self;
    typedef typename BaseInserter::size_type               size_type;
    typedef update_proxy<BaseInserter, size_type>          proxy_type;

    /// Constructor with matrix \p A, the mapping, and the slot size
    mapped_inserter(BaseInserter& A, Mapper& map)
      : ins(A), map(map) {}

  public:
    /// To be used in ins(r, c) << value;
    proxy_type operator[] (size_type row)        
    {	return proxy_type(ins, map.row(row));    }

    /// To be used in ins(r, c) << value;
    proxy_type operator() (size_type row)
    {	return proxy_type(ins, map.row(row));    }

    // update, modify and operator<< are used from BaseInserter

private:
    BaseInserter ins;
    Mapper&     map;
};

} // namespace vector
} // namespace mtl
#endif //MTL_VECTOR_MAPPED_INSERTER_INCLUDE
