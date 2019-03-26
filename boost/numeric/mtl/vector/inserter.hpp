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

#ifndef MTL_VECTOR_INSERTER_INCLUDE
#define MTL_VECTOR_INSERTER_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/operation/update.hpp>

namespace mtl { namespace vec {

/// Class for generic insertion into vectors
/** \param Vector is the type of the vector in which will be inserted
    \param Updater is a functor that determine how existed entries are updated
    \b Remark: For regular vectors, one can set the entries directly without inserter.
    An inserter is needed when vectors are distributed.
    Reversely, setting vectors with an inserter works for both distributed and non-distributed
    vector types and is therefore more generic. 
    To understand how to use inserters it is best to read the page \ref matrix_insertion.
    This applies all to vector insertion except that only one index is used.
**/
template <typename Vector, typename Updater>
struct inserter 
{
    typedef inserter                                 self;
    typedef typename Collection<Vector>::size_type   size_type;
    typedef typename Collection<Vector>::value_type  value_type;
    typedef update_proxy<self, size_type>            proxy_type;

    explicit inserter(Vector& vector) : vector(vector) {}

    proxy_type operator() (size_type n) { return proxy_type(*this, n); }
    proxy_type operator[] (size_type n) { return proxy_type(*this, n); }

    template <typename Modifier>
    void modify(size_type n, value_type value) { Modifier()(vector[n], value); }

    void update(size_type n, value_type value) { modify<Updater>(n, value); }

protected:
    Vector&    vector;
};


// Not to confuse with the equally-named class in operations (which should be in matrix anyway)
template <typename Inserter, typename Size>
struct update_proxy
{
    typedef update_proxy          self;
    typedef typename Inserter::value_type  value_type;

    explicit update_proxy(Inserter& ins, Size n) : ins(ins), n(n) {}
    
    template <typename Value>
    self& operator<< (Value const& val)
    {
	ins.update(n, val);
	return *this;
    }

    template <typename Value>
    self& operator= (Value const& val)
    {
	ins.template modify<update_store<value_type> > (n, val);
	return *this;
    }

    template <typename Value>
    self& operator+= (Value const& val)
    {
	ins.template modify<update_plus<value_type> > (n, val);
	return *this;
    }

  private:
    Inserter&  ins;
    Size       n;
};

}} // namespace mtl::vector

#endif // MTL_VECTOR_INSERTER_INCLUDE
