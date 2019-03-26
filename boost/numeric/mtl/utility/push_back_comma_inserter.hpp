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

#ifndef MTL_PUSH_BACK_COMMA_INSERTER_INCLUDE
#define MTL_PUSH_BACK_COMMA_INSERTER_INCLUDE

namespace mtl {

    /// Helper class to inserter with push_back using comma separation
    template <typename T>
    class push_back_comma_inserter
    {
	typedef push_back_comma_inserter self;
      public:
	/// Constructor takes a mutable reference of the object inserted into
	push_back_comma_inserter(T& ref) : ref(ref) {}

	/// Overloaded comma operator performs push_back
	template <typename Source>
	self& operator, (Source val)
	{ ref.push_back(val); return *this; }

      private:
	T& ref;
    };

} // namespace mtl

#endif // MTL_PUSH_BACK_COMMA_INSERTER_INCLUDE
