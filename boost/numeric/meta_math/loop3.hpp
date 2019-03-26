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

#ifndef META_MATH_LOOP3_INCLUDE
#define META_MATH_LOOP3_INCLUDE

// See below for example

namespace meta_math {

template <std::size_t Index0, std::size_t Max0, std::size_t Index1, std::size_t Max1,
	  std::size_t Index2, std::size_t Max2>
struct loop3
{
    static std::size_t const index0= Index0 - 1, next_index0= Index0,
	                     index1= Index1 - 1, next_index1= Index1,
            	             index2= Index2 - 1, next_index2= Index2 + 1;
};


template <std::size_t Index0, std::size_t Max0, std::size_t Index1, std::size_t Max1, 
	  std::size_t Max2>
struct loop3<Index0, Max0, Index1, Max1, Max2, Max2>
{
    static std::size_t const index0= Index0 - 1, next_index0= Index0,
	                     index1= Index1 - 1, next_index1= Index1 + 1,
            	             index2= Max2 - 1, next_index2= 1;
};


template <std::size_t Index0, std::size_t Max0, std::size_t Max1, std::size_t Max2>
struct loop3<Index0, Max0, Max1, Max1, Max2, Max2>
{
    static std::size_t const index0= Index0 - 1, next_index0= Index0 + 1,
	                     index1= Max1 - 1, next_index1= 1,
            	             index2= Max2 - 1, next_index2= 1;
};


template <std::size_t Max0, std::size_t Max1, std::size_t Max2>
struct loop3<Max0, Max0, Max1, Max1, Max2, Max2>
{
    static std::size_t const index0= Max0 - 1,
	                     index1= Max1 - 1,
            	             index2= Max2 - 1;
};




#if 0

// ============================
// Use the meta loop like this:
// ============================


template <std::size_t Index0, std::size_t Max0, std::size_t Index1, std::size_t Max1,
	  std::size_t Index2, std::size_t Max2>
struct loop3_trace : public loop3<Index0, Max0, Index1, Max1, Index2, Max2>
{
    typedef loop3<Index0, Max0, Index1, Max1, Index2, Max2> base;
    typedef loop3_trace<base::next_index0, Max0, base::next_index1, Max1, base::next_index2, Max2> next_t;

    void operator() ()
    {
	std::cout << this->index0 << " : " << this->index1 << " : " << this->index2 << "\n";
	next_t() ();
    }  
};


template <std::size_t Max0, std::size_t Max1, std::size_t Max2>
struct loop3_trace<Max0, Max0, Max1, Max1, Max2, Max2>
    : public loop3<Max0, Max0, Max1, Max1, Max2, Max2>
{
    void operator() ()
    {
	std::cout << this->index0 << " : " << this->index1 << " : " << this->index2 << "\n";
    }  
};

#endif

} // namespace meta_math

#endif // META_MATH_LOOP3_INCLUDE
