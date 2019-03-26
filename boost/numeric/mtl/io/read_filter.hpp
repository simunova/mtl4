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

#ifndef MTL_IO_READ_FILTER_INCLUDE
#define MTL_IO_READ_FILTER_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>

namespace mtl { namespace io {

/// Utility to filter entries in the read process
/** Depending on type only certain entries are considered for insertion.
    Particularly interesting for distributed collections (inserters). **/
template <typename Inserter>
class read_filter
{
  public:
    explicit read_filter(const Inserter& inserter) : inserter(inserter) {}
    
    /// Default for vectors is to consider every entry
    bool operator()(std::size_t) const { return true; }

    /// Default for matrices is to consider every entry
    bool operator()(std::size_t, std::size_t) const { return true; }

  private:
    const Inserter& inserter;
};


}} // namespace mtl::io

#endif // MTL_IO_READ_FILTER_INCLUDE
