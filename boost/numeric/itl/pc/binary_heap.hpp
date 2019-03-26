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
// Algorithm inspired by Nick Vannieuwenhoven, rewritten by Cornelius Steinhardt
//
//   Created on: Jan 10, 2010
//       Author: heazk


#ifndef MTL_BINARY_HEAP_INCLUDE
#define MTL_BINARY_HEAP_INCLUDE

namespace utils {

/**
 * An intrusive binary min-heap.
 *
 * KeyType: 		The type of the keys. Assumed to be a signed integer type.
 * KeyComparator:	A functor to compare values of the key type. Should be a
 * 					total order.
 */
template<
	class DirectAccessIterator,
	class KeyType,
	class KeyComparator,
	class ValueType,
	class GetKey,
	class GetParent,
	class GetLeft,
	class GetRight
>
class binary_heap {


/*******************************************************************************
 * Type Definitions
 ******************************************************************************/


public:
	/**
	 * The value type.
	 */
	typedef ValueType value_type;

	/**
	 * The key type.
	 */
	typedef KeyType key_type;

	/**
	 * The type of this heap.
	 */
	typedef binary_heap<
		DirectAccessIterator,
		KeyType,
		KeyComparator,
		ValueType,
		GetKey,
		GetParent,
		GetLeft,
		GetRight
	> heap_type;


/*******************************************************************************
 * Constructors and Destructors
 ******************************************************************************/


public:
	binary_heap(
		DirectAccessIterator values_begin,
		DirectAccessIterator values_end,
		KeyComparator compare_keys,
		GetKey get_key,
		GetParent get_parent,
		GetLeft get_left_child,
		GetRight get_right_child
	) :
		m_root(values_begin),
		m_compare_keys(compare_keys),
		m_key(get_key),
		m_parent(get_parent),
		m_left_child(get_left_child),
		m_right_child(get_right_child)
	{
		build_heap(values_begin, values_end);
	}

public:
	binary_heap(
		KeyComparator compare_keys,
		GetKey get_key,
		GetParent get_parent,
		GetLeft get_left_child,
		GetRight get_right_child
	) :
		m_root(0),
		m_last(0),
		m_compare_keys(compare_keys),
		m_key(get_key),
		m_parent(get_parent),
		m_left_child(get_left_child),
		m_right_child(get_right_child)
	{
	}

	/**
	 * Default destructor.
	 */
public:
	~binary_heap() {
		// Nothing here.
	}

	/**
	 * Copying of heaps is disallowed.
	 */
private:
	binary_heap(const heap_type&);
	void operator=(const heap_type&);


/*******************************************************************************
 * Heap Construction and Maintaining
 ******************************************************************************/


	/**
	 * Builds the min-heap.
	 */
private:
	void build_heap(
		DirectAccessIterator values_begin,
		DirectAccessIterator values_end
	) {
		std::cout << "building heap ..." << std::endl;
		// Construct heap connections.
		const int size = values_end - values_begin;
		m_parent(values_begin) = 0;
		m_left_child(values_begin) = 	1 < size ? values_begin+1 : 0;
		m_right_child(values_begin) = 	2 < size ? values_begin+2 : 0;
		for(int i = 1; i < size; ++i) {
			m_parent(values_begin+i) = values_begin+(((i+1)/2) - 1);
			m_left_child(values_begin+i) =
				2*(i+1)-1 < size ? values_begin+(2*(i+1)-1) : 0;
			m_right_child(values_begin+i) =
				2*(i+1) < size ? values_begin+(2*(i+1)) : 0;
		}

		m_root = size > 0 ? values_begin : 0;
		m_last = size > 0 ? values_end-1 : 0;

		build_heap_rec(m_root);

	}

	/**
	 * Applies the min-heapify procedure at each node in postfix style.
	 */
private:
	void build_heap_rec(value_type node) {
		if(node) {
//			dump(node);
			build_heap_rec(m_left_child(node));
			build_heap_rec(m_right_child(node));
			min_heapify(node);
		}
	}

	/**
	 * The min-heapify operation to restore the min heap property at a given
	 * node.
	 */
private:
	void min_heapify(value_type node) {

//		std::cout << "heapify: " << m_key(node) << std::endl;

		value_type seek = node;
		value_type smallest = seek;
		do {
			value_type left_child = m_left_child(seek);
			value_type right_child = m_right_child(seek);
			smallest = seek;

			if(
					left_child &&
					m_compare_keys(left_child, seek)
			) {
				smallest = left_child;
			}
			if(
					right_child &&
					m_compare_keys(right_child, smallest)
			) {
				smallest = right_child;
			}
			if(smallest != seek) {
				exchange(smallest, seek);
			}
		} while(smallest != seek);

	}

	/**
	 * Exchanges two nodes.
	 */
private:
	inline void exchange(value_type fst, value_type snd) {
		if(!fst && !snd) {
			return;
		}
		if(!fst || !snd) {
			assert(false);
		}
		if(fst == snd) {
			return;
		}

		// Update root and last pointers.
		if(fst == m_root) {
			m_root = snd;
		} else if(snd == m_root) {
			m_root = fst;
		}
		if(fst == m_last) {
			m_last = snd;
		} else if(snd == m_last) {
			m_last = fst;
		}

//		std::cout << "fst pre" << std::endl;
//		dump(fst);
//		std::cout << "snd pre" << std::endl;
//		dump(snd);

		// Canonize the situation when fst and snd are directly connected such
		// that fst is always the higher node.
		if( m_parent(fst) == snd ) {
			value_type tmp = fst;
			fst = snd;
			snd = tmp;
		}

		// Set the new parent pointers for the children.
		value_type left_child = m_left_child(fst);
		value_type right_child = m_right_child(fst);
		if(left_child) {
			m_parent(left_child) = snd;
		}
		if(right_child) {
			m_parent(right_child) = snd;
		}
		left_child = m_left_child(snd);
		right_child = m_right_child(snd);
		if(left_child) {
			m_parent(left_child) = fst;
		}
		if(right_child) {
			m_parent(right_child) = fst;
		}

		// Set the new children pointers.
		const value_type fst_left_child = m_left_child(fst);
		const value_type fst_right_child = m_right_child(fst);
		m_left_child(fst) = m_left_child(snd);
		m_right_child(fst) = m_right_child(snd);
		m_left_child(snd) = fst_left_child;
		m_right_child(snd) = fst_right_child;

		// Set the new children pointers of the parents.
		if(m_parent(fst)) {
			if( m_left_child(m_parent(fst)) == fst ) {
				m_left_child(m_parent(fst)) = snd;
			} else {
				m_right_child(m_parent(fst)) = snd;
			}
		}
		if(m_parent(snd)) {
			if( m_left_child(m_parent(snd)) == snd ) {
				m_left_child(m_parent(snd)) = fst;
			} else {
				m_right_child(m_parent(snd)) = fst;
			}
		}

		// Set the new parent pointers.
		const value_type fst_parent = m_parent(fst);
		m_parent(fst) = m_parent(snd);
		m_parent(snd) = fst_parent;

//		std::cout << "fst post" << std::endl;
//		dump(snd);
//		std::cout << "snd post" << std::endl;
//		dump(fst);

//		if( m_parent(fst) == snd || m_parent(snd) == fst) {
//			// Canonize the situation to this whereby fst is the higher node.
//			if( m_parent(fst) == snd ) {
//				value_type tmp = fst;
//				fst = snd;
//				snd = tmp;
//			}
//
//			// Update root and last pointers.
//			if(fst == m_root) {
//				m_root = snd;
//			}
//			if(snd == m_last) {
//				m_last = fst;
//			}
//
//			// Get a hold on every value we need.
//			value_type fst_parent = m_parent(fst);
//			value_type fst_left_child = m_left_child(fst);
//			value_type fst_right_child = m_right_child(fst);
//			value_type snd_left_child = m_left_child(snd);
//			value_type snd_right_child = m_right_child(snd);
//
//			// Set up fst's new pointers.
//			m_left_child(fst) = snd_left_child;
//			m_right_child(fst) = snd_right_child;
//			m_parent(snd_left_child) = fst;
//			m_parent(snd_right_child) = fst;
//			m_parent(fst) = snd;
//
//			// Set up snd's new pointers.
//			m_parent(snd) = fst_parent;
//			if(fst_parent) {
//				if( m_left_child(fst_parent) == fst ) {
//					m_left_child(fst_parent) = snd;
//				} else {
//					m_right_child(fst_parent) = snd;
//				}
//			}
//
//			if(fst_left_child == snd) {
//				m_parent(fst_right_child) = snd;
//				m_right_child(snd) = fst_right_child;
//				m_left_child(snd) = fst;
//			} else {
//				m_parent(fst_left_child) = snd;
//				m_left_child(snd) = fst_left_child;
//				m_right_child(snd) = fst;
//			}
//		}
	}


/*******************************************************************************
 * Heap Operations
 ******************************************************************************/

	/**
	 * Checks whether the heap is empty.
	 */
public:
	bool empty() {
		return m_root == 0;
	}

	/**
	 * Returns a reference to the top of the heap.
	 */
public:
	value_type& top() {
		return m_root;
	}

	/**
	 * Pops the top of the heap.
	 */
public:
	void pop() {
		assert(top());
		assert(m_last);

//		dump(m_root);
//		dump(m_last);

		value_type new_last = m_last;
		if(m_root != m_last) {

//			std::cout << "gonna	 exchange" << std::endl;
			exchange(m_root, m_last);

			assert(m_root == new_last);

			new_last = m_last;
//			std::cout << "donna	 exchange" << std::endl;

//			dump(m_root);
//			dump(m_last);

			// Determine the new "last" node of the tree.
			// The algorithm also works when the last node at the deepest level
			// is removed. The algorithm will then continue from the rightmost
			// node at the level one lower.
			while(
					m_parent(new_last) &&
					m_left_child(m_parent(new_last)) == new_last
			) {
				new_last = m_parent(new_last);
//				std::cout << "moving up ..." << std::endl;
			}
//			std::cout << "stopped moving up ..." << std::endl;
			if(m_parent(new_last)) {
				new_last = m_left_child(m_parent(new_last));
//				std::cout << "moved left ..." << std::endl;
			}
			while(m_right_child(new_last)) {
				new_last = m_right_child(new_last);
//				std::cout << "moving right ..." << std::endl;
			}
//			std::cout << "stopped moving right ..." << std::endl;
		}

		if( m_parent(m_last) ) {
			if( m_left_child(m_parent(m_last)) == m_last ) {
				m_left_child(m_parent(m_last)) = 0;
			} else {
				m_right_child(m_parent(m_last)) = 0;
			}
		}
		m_parent(m_last) = 0;
		m_left_child(m_last) = 0;
		m_right_child(m_last) = 0;
//		m_key(m_last) = -1;

		if(m_root != m_last) {
			m_last = new_last;
			min_heapify(m_root);
		} else {
			m_root = m_last = 0;
		}
	}

	/**
	 * Updates the key of a given node.
	 */
public:
	void update_key(value_type node, key_type new_key) {
		m_key(node) = new_key;

		value_type seek = node;
		while(
				m_parent(seek) &&
				m_compare_keys(seek, m_parent(seek))
		) {
			exchange(seek, m_parent(seek));
		}
	}

	/**
	 * Inserts a new value into the binary heap with the given key.
	 */
public:
	void insert(value_type node, key_type key) {

		// Special case handling an empty heap.
		if(m_last == m_root && m_last == 0) {
			m_last = m_root = node;
			m_parent(node) = 0;
			m_left_child(node) 	= 0;
			m_right_child(node) = 0;
			m_key(node) = key;
			return;
		}

		// Search the new "last" position to insert the new node at.
		value_type new_last = m_last;

		// While the current node is right w.r.t. its parent, move up the tree.
		while(
				m_parent(new_last) &&
				m_right_child(m_parent(new_last)) == new_last
		) {
			new_last = m_parent(new_last);
		}
		// If there is a right sibling, go to it and descent to the left.
		if(m_parent(new_last) && m_right_child(m_parent(new_last))) {
			new_last = m_right_child(m_parent(new_last));
			while(m_left_child(new_last)) {
				new_last = m_left_child(new_last);
			}
		}
		// If there is a parent at this point, the right sibling does not exist,
		// and the given node should be placed there.
		if(m_parent(new_last)) {
			// Add the node to the tree.
			if(m_right_child(m_parent(new_last)) == 0) {
				new_last = m_parent(new_last);
				m_right_child(new_last) = node;
			} else {
				m_left_child(new_last) = node;
			}
		} else {
			// We came from the very last node at the current level and hence
			// must begin a new level. Descent to the left.
			while( m_left_child(new_last) ) {
				new_last = m_left_child(new_last);
			}

			// Add the node to the tree.
			m_left_child(new_last) = node;
		}

		m_parent(node) = new_last;
		m_left_child(node) 	= 0;
		m_right_child(node) = 0;
		m_key(node) = key;
		m_last = node;

		// Update the key and force the heap to rebalance itself.
		update_key(node, key);
	}

	/**
	 * Removes the node. The operation can safely be applied to an element no
	 * longer in the heap.
	 */
public:
	void remove(value_type node) {
		if( m_parent(node) == 0 && m_root != node) {
			return;
		}

		// Move the node up the tree.
		value_type seek = node;
		while(m_parent(seek)) {
			exchange(seek, m_parent(seek));
		}
		// Remove it.
		pop();
	}


public:
	void dump(value_type node) {
		std::cout << "seq: " << node->get_sequence_number() << std::endl;
		std::cout << "key: " << m_key(node) << std::endl;
		std::cout << "par: " << (m_parent(node) ? (m_parent(node))->get_sequence_number() : -1) << std::endl;
		std::cout << "lef: " << (m_left_child(node) ? (m_left_child(node))->get_sequence_number() : -1) << std::endl;
		std::cout << "rig: " << (m_right_child(node) ? (m_right_child(node))->get_sequence_number() : -1) << std::endl;

	}

/*******************************************************************************
 * Data Members
 ******************************************************************************/


public:
	/**
	 * The root of the binary heap.
	 */
	ValueType m_root;

	/**
	 * The "last" node of the binary heap. The last node is the node that can
	 * safely be removed such that the resulting tree is still a complete binary
	 * tree, with the additional constraint that the nodes in the last level are
	 * all in the "left" part of the tree.
	 */
	ValueType m_last;

	/**
	 * A functor to compare keys.
	 */
	KeyComparator m_compare_keys;

	/**
	 * A functor to obtain the key.
	 */
	GetKey m_key;

	/**
	 * A functor to obtain the parent.
	 */
	GetParent m_parent;

	/**
	 * A functor to obtain the left child.
	 */
	GetLeft m_left_child;

	/**
	 * A functor to obtain the right child.
	 */
	GetRight m_right_child;

};

}

#endif // MTL_BINARY_HEAP_INCLUDE
