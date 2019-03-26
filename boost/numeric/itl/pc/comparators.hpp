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
//
// Algorithm inspired by Nick Vannieuwenhoven, written by Cornelius Steinhardt



/*
 *
 *  Created on: Oct 14, 2010
 *      Author: heazk
 */


#ifndef MTL_COMPARATORS_INCLUDE
#define MTL_COMPARATORS_INCLUDE

#define POINTER long

namespace compare {

/**
 * Compares two given pointers based on the address to which they point.
 */
template<class Type>
struct address_compare {
	bool operator()(const Type *const a, const Type *const b) const {
		return reinterpret_cast<POINTER>(a) < reinterpret_cast<POINTER>(b);
	}
};

/**
 * Compares two given pointers based on the address to which they point.
 */
template<class Type>
struct address_compare_equal {
	bool operator()(const Type *const a, const Type *const b) const {
		return reinterpret_cast<POINTER>(a) == reinterpret_cast<POINTER>(b);
	}
};

/**
 * Hashes a given pointer based on its address.
 */
template<class Type>
struct address_hasher {
	POINTER operator()(const Type *const a) const {
		return reinterpret_cast<POINTER>(a);
	}
};

} // end namespace compare


#endif // MTL_COMPARATORS_INCLUDE
