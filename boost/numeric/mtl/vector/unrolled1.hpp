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

#ifndef MTL_VECTOR_UNROLLED1_INCLUDE
#define MTL_VECTOR_UNROLLED1_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/crtp_base_vector.hpp>

namespace mtl { namespace vec {

/// Helper class for unrolling loops \sa \ref mtl::unroll
template <unsigned BSize, typename Vector>
class unrolled1
  : public crtp_base_vector< unrolled1<BSize, Vector>, typename Collection<Vector>::value_type, std::size_t >
{
    typedef unrolled1                               self;
  public:
    typedef typename Collection<Vector>::value_type value_type;
    typedef typename Collection<Vector>::size_type  size_type;
    typedef value_type&                             reference ;
    typedef crtp_vector_assign< self, value_type, std::size_t > assign_base; // base of crtp_base_vector

    unrolled1(Vector& ref) : ref(ref) {}
	    
    reference operator()(size_type i) { return ref(i); }
    reference operator[](size_type i) { return (*this)(i); }
    // const versions shouldn't be needed because it is supposed to be a lvalue

    //friend inline size_type size(const self& v)  { return size(v.ref); }
    friend inline size_type num_rows(const self& v)  { return num_rows(v.ref); }
    friend inline size_type num_cols(const self& v)  { return num_cols(v.ref); }

    void change_dim(size_type d) { ref.change_dim(d); }

    using assign_base::operator=;

#if !defined(_MSC_VER) || _MSC_VER != 1400 // Bug in MSVC 2005
    template <unsigned BBSize, typename VVector>
    friend inline std::size_t size(const unrolled1<BBSize, VVector>& v);
  private:
#endif
    Vector&    ref;
};

template <unsigned BSize, typename Vector>
inline std::size_t size(const unrolled1<BSize, Vector>& v)
{ return size(v.ref); }

}} // namespace mtl::vector

#endif // MTL_VECTOR_UNROLLED1_INCLUDE
