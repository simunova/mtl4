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

#ifndef MTL_BASE_CURSOR_INCLUDE
#define MTL_BASE_CURSOR_INCLUDE

namespace mtl { namespace detail {

// base class for different cursors, works with pointers and integers
template <class Key> class base_cursor 
{
 public:
    typedef Key          key_type;
    typedef base_cursor  self;

    base_cursor () {} 
    base_cursor (key_type key) : key(key) {}

    key_type operator*() const 
    { 
      return key; 
    }

    key_type value() const 
    { 
      return key; 
    }

    self& operator++ () 
    { 
      ++key; return *this; 
    }
    self operator++ (int) 
    { 
      self tmp = *this; 
      ++key; 
      return tmp; 
    }
    self& operator-- () 
    { 
      --key; 
      return *this; 
    }
    self operator-- (int) 
    { 
      self tmp = *this; 
      --key; 
      return tmp; 
    }
    template <typename T>
    self& operator+=(T n) 
    { 
      key += key_type(n); 
      return *this; 
    }
  
    template <typename T>
    self operator+(T n) const
    {
	self tmp = *this;
	tmp+= n;
	return tmp;
    }
    
    template <typename T>
    self& operator-=(T n)
    { 
	key -= key_type(n); 
	return *this; 
    }

    template <typename T>
    self operator-(T n) const
    {
	self tmp = *this;
	tmp -= n;
	return tmp;
    }

    int operator-(const self& cc) const
    {
	return this->key - cc.key;
    }

    bool operator==(const self& cc) const 
    {
      return key == cc.key; 
    }

    bool operator!=(const self& cc) const 
    {
      return !(*this == cc); 
    }
  
  


    key_type key;
}; // base_cursor



}} // namespace mtl::detail 

#endif // MTL_BASE_CURSOR_INCLUDE 


