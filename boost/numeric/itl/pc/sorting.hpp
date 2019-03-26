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
 *  Created on: Nov 3, 2009
 *      Author: heazk
 */

#ifndef MTL_SORTING_INCLUDE
#define MTL_SORTING_INCLUDE

#include <string.h>
#include <stdlib.h>

#define BITS_PER_BYTE 8

////////////////////////////////////////////////////////////////////////////////
// RADIX-SORT
////////////////////////////////////////////////////////////////////////////////

template<
	class Type
>
void radix_sort(
	Type* values,
	const int size,
	Type* tmp = 0
) {
	if(size < 2) {
		return;
	}

	Type *const orig_ptr = values;
	if(tmp == 0) {
		tmp = new Type[size];
	}

#if VERBOSE_MODE > 10
	std::cout << "unsorted: ";
	for(int i = 0; i < size; ++i) {
		std::cout << values[i] << ", ";
	}
	std::cout << std::endl;
#endif

	// Determine the maximum number of bits needed.
	int max = INT_MIN;
	for(int i = 0; i < size; ++i) {
		max = (values[i] <= max) ? max : values[i];
	}
	int maxBits = 1;
	for(; max; max >>= 1, ++maxBits){};
	for(max = 0; max < maxBits; max += BITS_PER_BYTE){};
	maxBits = max;

	const int groupBits = BITS_PER_BYTE;
	const int intBits   = maxBits;
	const int groupSize = 1 << groupBits;

	const int nbIterations = intBits / groupBits;
	const int mask = groupSize-1;

	// Counting and prefix arrays.
	int* count = new int[groupSize];

	for (
			int shift = 0, it = 0;
			it < nbIterations;
			++it, shift += groupBits
	) {
		// Reset count array.
		memset(count, 0, groupSize*sizeof(int));

		// Counting elements in the it-th group.
		for (int i = 0; i < size; ++i) {
			count[ (values[i] >> shift) & mask ]++;
		}

		// Calculate prefixes.
		int currIncrement = -1;
		int prevCount = count[0];
		count[0] = 0;
		for (int i = 1; i < groupSize; ++i) {
			currIncrement = prevCount;
			prevCount = count[i];
			count[i] = count[i-1] + currIncrement;
		}

		// Sort elements in the it-th group.
		for (int i = 0; i < size; ++i) {
			int& index = count[( values[i] >> shift) & mask];
			tmp[ index ] = values[i];
			++index;
		}

		// Swap pointers.
		Type* tmp_values = values;
		values = tmp;
		tmp = tmp_values;
	}

	// Make sure the sorted keys and values reside in the input array.
	if( values != orig_ptr ) {
		memcpy( orig_ptr, values, sizeof(Type)*size );
		tmp = values;
	}

#if VERBOSE_MODE > 10
	std::cout << "sorted: ";
	for(int i = 0; i < size; ++i) {
		std::cout << values[i] << ", ";
	}
	std::cout << std::endl;
#endif

	if(count) {
		delete[] count;
	}
	if(tmp) {
		delete[] tmp;
	}
}

template<
	class Key, class Val
>
void radix_sort(
		Key* keys, Val* vals,
		const int size,
		Key* tmpKeys = 0, Val* tmpVals = 0
) {
	if(size < 2) {
		return;
	}

	// The pointers to the original data.
	Key *const origKey = keys;
	Val *const origVal = vals;

	// Helper arrays.
	bool ownKeys = false;
	bool ownVals = false;
	if(!tmpKeys) {
		tmpKeys = new Key[size];
		ownKeys = true;
	}
	if(!tmpVals) {
		tmpVals = new Val[size];
		ownVals = true;
	}

	// Determine the maximum number of bits needed.
	int max = -(INT_MAX-1);
	for(int i = 0; i < size; ++i) {
		max = keys[i] <= max ? max : keys[i];
	}
	int maxBits = 1;
	for(; max; max >>= 1, ++maxBits){};
	for(max = 0; max < maxBits; max += BITS_PER_BYTE){};
	maxBits = max;

	const int groupBits = BITS_PER_BYTE;
	const int intBits   = maxBits;
	const int groupSize = 1 << groupBits;

	const int nbIterations = intBits / groupBits;
	const int mask = groupSize-1;

	// Counting and prefix arrays.
	int* count = new int[groupSize];
	// Allocation of memory could fail!

	for (
			int shift = 0, it = 0;
			it < nbIterations;
			++it, shift += groupBits
	) {
		// Reset count array.
		for(int i = 0; i < groupSize; ++i) {
			count[i] = 0;
		}

		// Counting elements in the it-th group.
		for (int i = 0; i < size; ++i) {
			count[ (keys[i] >> shift) & mask ]++;
		}

		// Calculate prefixes.
		int currIncrement = -1;
		int prevCount = count[0];
		count[0] = 0;
		for (int i = 1; i < groupSize; ++i) {
			currIncrement = prevCount;
			prevCount = count[i];
			count[i] = count[i-1] + currIncrement;
		}

		// Sort elements in the it-th group.
		for (int i = 0; i < size; ++i) {
			int& index = count[( keys[i] >> shift) & mask];
			tmpKeys[ index ] = keys[i];
			tmpVals[ index ] = vals[i];
			++index;
		}

		// Swap pointers.
		Key* tmpKeyPtr = tmpKeys;
		Val* tmpValuePtr = tmpVals;
		tmpKeys = keys;
		tmpVals = vals;
		keys = tmpKeyPtr;
		vals = tmpValuePtr;
	}

	// Make sure the sorted keys and values reside in the input array.
	if(keys != origKey) {
		memcpy( origKey, keys, sizeof(Key)*size );
		tmpKeys = keys;
	}
	if(vals != origVal) {
		memcpy( origVal, vals, sizeof(Val)*size );
		tmpVals = vals;
	}

	if(count) 				delete[] count;
	if(ownKeys && tmpKeys) 	delete[] tmpKeys;
	if(ownVals && tmpVals)	delete[] tmpVals;
	tmpKeys = 0;
	tmpVals = 0;
	count = 0;
}

////////////////////////////////////////////////////////////////////////////////
// QUICK-SORT
////////////////////////////////////////////////////////////////////////////////

template<typename Key, typename Val>
inline void swap(Key& key0, Val& val0, Key& key1, Val& val1) {
	Key tmp_key = key0;
	Val tmp_val = val0;
	key0 = key1;
	val0 = val1;
	key1 = tmp_key;
	val1 = tmp_val;
}

template<typename Key, typename Value>
inline void sort_along(Key* keys, Value* values, int size) {


	if(size < 2) {
		return;
	}
	if( size == 2 ) {
		if(keys[0] > keys[1]) {
			swap(
				keys[0], values[0],
				keys[1], values[1]
			);
		}
		return;
	}

	int piv = std::rand() % size;


	swap(
		keys[0],   values[0],
		keys[piv], values[piv]
	);
	int begin = 1;
	int end = size-1;
	while(begin < end) {
		for(; (begin < end) && (keys[begin] <= keys[0]); ++begin){};
		for(; (begin <= end) && (keys[end] > keys[0]); --end){};

		if(begin < end) {
			swap(
				keys[begin], values[begin],
				keys[end],   values[end]
			);
		}
	}
	piv = end;
	swap(
		keys[0],   values[0],
		keys[piv], values[piv]
	);


	sort_along( keys, values, piv );
	sort_along( keys+(piv+1), values+(piv+1), size-piv-1);
}

#endif // MTL_SORTING_INCLUDE
