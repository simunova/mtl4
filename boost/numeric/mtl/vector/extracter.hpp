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
// See also license.short.txt in the distribution.

#ifndef MTL_VECTOR_EXTRACTER_INCLUDE
#define MTL_VECTOR_EXTRACTER_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>

namespace mtl {	namespace vec {

template< typename Vector >
struct indexbuffer 
{
    typedef typename Collection<Vector>::value_type value_type;
    typedef typename Collection<Vector>::size_type size_type;
    typedef indexbuffer self;

    indexbuffer(const Vector& v) : vec(v) {}

    inline self& operator()(const size_type p, value_type& dest) 
    { 
	dest = vec[p];
	return *this; 
    }

    inline void update() {}
  private:
    const Vector& vec;
};

	
template< typename Vector >
class extracter 
{			
  public:
    typedef typename Collection<Vector>::value_type value_type;
    typedef typename Collection<Vector>::size_type  size_type;
    typedef indexbuffer< Vector >                   buffer_type;

    extracter(const Vector& vec):
      destroyBuffer(true), buffer(new buffer_type(vec))
    {}
    extracter(buffer_type& buf):
      destroyBuffer(false), buffer(&buf) 
    {}
			
    struct bracket_proxy 
    {
	const size_type p;
	buffer_type& buffer;
	bracket_proxy(const size_type p, buffer_type& b) : p(p), buffer(b) {}

	inline buffer_type& operator>>(value_type& dest) { return buffer(p, dest); }
    };

    inline bracket_proxy operator[](const size_type p) { return bracket_proxy(p, *buffer); }

    ~extracter() {
	buffer->update();
	if(destroyBuffer)
	    delete buffer;
    }
  private:
    const bool destroyBuffer;
    buffer_type* buffer;
};
	
}} //namespace mtl::vector
#endif
