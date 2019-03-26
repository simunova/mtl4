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

 
#ifndef MTL_ELEMENT_INCLUDE
#define MTL_ELEMENT_INCLUDE

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <functional>

#include <boost/unordered_set.hpp>

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/utility/make_copy_or_reference.hpp>
#include <boost/numeric/itl/pc/comparators.hpp>


namespace mtl { namespace mat {

///  A class representing an element with ValType equals the type of the numeric values
template <typename ValType>
class element 
{

public:
    /// The type of this element.
    typedef element<ValType> element_type;

    /// The value type of the matrix and rhs elements.
    typedef ValType value_type;

    /// The type of a set of neighbors.
    typedef std::vector<element_type*, std::allocator<element_type*> > neighbor_collection_type; // trouble with Clang 4.2 w/o allocator 

    /// An iterator over the neighbors of this element.
    typedef typename neighbor_collection_type::iterator neighbor_iterator;

    /// The type of an unordered set of neighbors.
    typedef typename boost::unordered_set<
	element_type*,
	compare::address_hasher<element_type>,
	compare::address_compare_equal<element_type>
       > neighbor_set_type;

    /// The type of the iterator over an unorderd set of neighbors.
    typedef typename neighbor_set_type::iterator neighbor_set_iterator_type;

    /// The type of matrix.
    typedef mtl::mat::dense2D<value_type> matrix_type;

    /// The type of the index vector.
    typedef mtl::dense_vector<int> index_type; 


/*******************************************************************************
 * Constructors
 ******************************************************************************/


    /**
     * Constructs the element using the memory specified by the two references.
     *
     * p_indices: 	a reference to the memory where the indices may be stored.
     * p_values:	a reference to the memory where the values may be stored.
     */
  public:
    element(int p_sequence_number, const index_type& p_indices, const matrix_type& p_values) 
      : m_indices(p_indices),	
	m_values(p_values),
	m_sequence_number(p_sequence_number)
    {}

    element() 
      :	m_sequence_number(-1) {}

    element(const element_type& other) 
      :	m_sequence_number(-1)
    {	*this = other;	}

    /// Deep copy the given element.
    void operator=(const element_type& other) 
    {
	m_sequence_number = other.m_sequence_number;
	m_neighbors = other.m_neighbors;

	m_indices = other.m_indices;
	m_values =  other.m_values;
    }

    /// Returns the unique identifier of this element.
    inline int get_id() const {	return m_sequence_number; }
    inline int& get_id()      {return m_sequence_number; }

    /// Returns the number of variables.
    inline std::size_t nb_vars() const { return size(m_indices); }

    /// Returns the number of values.
    inline std::size_t nb_values() const { return nb_vars()*nb_vars(); }

    /// Reference to the matrix of values.
    inline matrix_type& get_values() { return m_values; }

    /// Constant reference to the matrix of values.
    inline const matrix_type& get_values() const { return m_values; }

    /// Reference to the indices.
    inline const index_type& get_indices() const { return m_indices; }

    /// Mutable reference to the indices.
    inline index_type& get_indices() {	return m_indices; }

    /// The actual number of non-zero values.
    int nnz() const 
    {
	const value_type    zero= math::zero(value_type());
	int nbr_nz = 0;
	for (int r = 0; r < nb_vars(); ++r) 
	    for (int c = 0; c < nb_vars(); ++c) 
		nbr_nz +=  get_values()(r,c) != zero;
	return nbr_nz;
    }

    /// Reference to the set of neighbors.
    neighbor_collection_type& get_neighbors() { return m_neighbors;	}

    /// Reference to the set of neighbors.
    const neighbor_collection_type& get_neighbors() const { return m_neighbors; }

    /// Number of neighbors this element is connected to.
    int get_nb_neighbors() const {	return int(m_neighbors.size()); }

    /// Add new neighbors, max 6 at the time
    void add_neighbors(element* n1, element* n2= 0, element* n3= 0,
			element* n4= 0, element* n5= 0, element* n6= 0) 
    { 
	m_neighbors.push_back(n1);
	if (n2) {
	    m_neighbors.push_back(n2);
	    if (n3) m_neighbors.push_back(n3);
	    if (n4) m_neighbors.push_back(n4);
	    if (n5) m_neighbors.push_back(n5);
	    if (n6) m_neighbors.push_back(n6);
	}
    }


/*******************************************************************************
 * Useful Inspector Methods
 ******************************************************************************/

    /// The set of nodes that is incident to the element.
    boost::unordered_set<int> get_incident_nodes() const 
    {
	boost::unordered_set<int> nodes(2 * get_nb_neighbors());
	for(typename neighbor_collection_type::const_iterator neigh_it = m_neighbors.begin();
	    neigh_it != m_neighbors.end(); ++neigh_it) {
	    element_type& neigh = **neigh_it;
	    nodes.insert(neigh.get_indices().begin(), neigh.get_indices().end());
	}
	// Remove the nodes of the element.
	for (std::size_t i = 0; i < nb_vars(); ++i ) 
	    nodes.erase( get_indices()(i) );
	return nodes;
    }


    /// Get the set of level-k neighbors, for a given k. 
    neighbor_set_type get_level_neighbors(const int level = 1) 
    {
	neighbor_set_type result( get_nb_neighbors() * level );

	if (level > 0) {
	    result.insert( m_neighbors.begin(), m_neighbors.end() );
	    if (level > 1) {
		for(int i = 0; i < get_nb_neighbors(); ++i) {
		    neighbor_set_type neighs(m_neighbors[i]->get_level_neighbors(level-1));
		    result.insert( neighs.begin(), neighs.end() );
		}
		result.erase( this );
	    }
	}
	return result;
    }

/*******************************************************************************
 * Manipulation
 ******************************************************************************/
  public:
    /// Permutes the rows and the columns of the element coefficient matrix along 
    /// with the indices such that the latter are sorted in ascending order.
    void sort_indices() {
	if(size(m_indices) == 0) {
	    assert(size(m_values) == 0);
	    return;
	}

	bool sorted = true;
	for (std::size_t i = 0; i < nb_vars()-1; ++i) 
	    sorted &= (get_indices()(i) < get_indices()(i+1));
	if (sorted) 
	    return;	

	index_type  orig_index( get_indices() );
	matrix_type orig_matrix( get_values() );

	std::sort(
		  &(get_indices()(0)),
		  &(get_indices()(0))+nb_vars()
		  );

	index_type orig_offset( nb_vars() );
	orig_offset = -1;
	for(std::size_t i = 0; i < nb_vars(); ++i) {
	    int seek_idx = get_indices()(i);
	    std::size_t j = 0;
	    for(; j < nb_vars() 
		    && orig_index(j) != seek_idx; ++j){};
	    orig_offset(i) = int(j);
	}

	matrix_type& values = get_values();
	for(std::size_t r = 0; r < nb_vars(); ++r) {
	    for(std::size_t c = 0; c < nb_vars(); ++c) {
		values(r,c) = orig_matrix( orig_offset(r), orig_offset(c) );
	    }
	}

#ifndef NDEBUG
	sorted = true;
	for(std::size_t i = 0; i < nb_vars()-1; ++i) {
	    sorted &= (get_indices()(i) < get_indices()(i+1));
	}
	assert(sorted);
#endif
    }


  public:
    /// Removes the given set of nodes from the element
    template< class Vector >
    void remove_nodes(const Vector& nodes, element_type& el) {
	if(size(m_indices) == 0) {
	    assert(size(m_values) == 0);
	    return;
	}
	if(nb_vars() == 0) {
	    return;
	}

#ifndef NDEBUG
	bool sorted = true;
	for(unsigned int i = 1; i < mtl::size(nodes); ++i) {
	    sorted &= ( nodes[i-1] < nodes[i] );
	}
	assert(sorted);
#endif

	const std::size_t nb_nodes = mtl::size(nodes);

	// Count number of remaining variables.
	long new_nb_nodes = long(nb_vars());
	{
	    std::size_t i = 0, j = 0;
	    while( i < nb_vars() && j < nb_nodes ) {
		const int diff = get_indices()(i) - nodes[j];
		if( diff < 0 ) {
		    ++i;
		} else if( diff > 0 ) {
		    ++j;
		} else {
		    --new_nb_nodes;
		    ++i;
		    ++j;
		}
	    }
	}
	assert(new_nb_nodes >= 0);

	// Construct new index array.
	index_type index;
	index_type local_index(new_nb_nodes);
	if (new_nb_nodes > 0) {
	    index.change_dim(new_nb_nodes);
	    std::size_t i = 0, j = 0, pos = 0;
	    while( i < nb_vars() && j < nb_nodes ) {
		const int diff = get_indices()(i) - nodes[j];
		if( diff < 0 ) {
		    assert( pos < std::size_t(new_nb_nodes) );
		    index[pos] = get_indices()(i);
		    local_index(pos) = int(i);
		    ++pos;
		    ++i;
		} else if( diff > 0 ) {
		    ++j;
		} else {
		    ++i;
		    ++j;
		}
	    }
	    while( i < nb_vars() ) {
		assert( pos < std::size_t(new_nb_nodes) );
		index[pos] = get_indices()(i);
		local_index(pos) = int(i);
		++pos;
		++i;
	    }
	} 

	matrix_type values;
	if(new_nb_nodes > 0) {
	    values.change_dim( new_nb_nodes, new_nb_nodes );
	    matrix_type tmp(get_values()), tmp2(new_nb_nodes, new_nb_nodes);
	    for(unsigned int i=0;i<size(local_index);i++){
		for(unsigned int j=0;j<size(local_index);j++){
		    tmp2[i][j]=tmp[local_index(i)][local_index(j)];
		}
	    }
	    values = tmp2;
	} else {
//	std::cout<< "ELSE\n";
//	    values = new matrix_type(0,0);
	}
	// Update the neighborhood.
	std::set<int, std::less<int>, std::allocator<int> > remove_neighs; // trouble with Clang 4.2 w/o allocator 
	for(
	    neighbor_iterator neigh_it = m_neighbors.begin();
	    neigh_it != m_neighbors.end();
	    ++neigh_it
	    ) {
	    element_type& neigh = **neigh_it;

	    // Search a matching index.
	    bool connected = false;
	    {
		std::size_t i = 0, j = 0;
		while(
		      i < std::size_t (new_nb_nodes) &&
		      j < neigh.nb_vars() &&
		      !connected
		      ) {
//		    const int diff = (*index)(i) - neigh.get_indices()(j);
		    const int diff = index[i] - neigh.get_indices()(j);
		    if(diff < 0) {
			++i;
		    } else if(diff > 0) {
			++j;
		    } else {
			connected = true;
		    }
		}
	    }

	    // If not found, then remove ourself from the neighbors and vice versa.
	    if(!connected) {
		neighbor_iterator pos = std::find(neigh.get_neighbors().begin(), neigh.get_neighbors().end(), this );
		if( (pos != neigh.get_neighbors().end()) && (&neigh != &el) ) {
		    neigh.get_neighbors().erase(pos);
		}
		remove_neighs.insert( neigh.get_id() );
	    }
	}

	// Remove the neighbors we're no longer connected to.
	for(std::set<int>::iterator it = remove_neighs.begin(); it != remove_neighs.end(); ++it) {
	    const int seek_seq_nbr = *it;
	    for (std::size_t j = 0; j < m_neighbors.size(); ++j) 
		if (m_neighbors[j] != 0 && m_neighbors[j]->get_id() == seek_seq_nbr) {
		    m_neighbors.erase( m_neighbors.begin()+j );
		    break;
		}
	}

	if(new_nb_nodes == 0) {
	    m_neighbors.clear();
	}

	m_indices.change_dim(0);
	m_values.change_dim(0, 0);
	m_indices = index;
	m_values = values;
    }

    /// Absorbs the values of the given matrix with the given index.
    template< class Matrix, class Vector >
    void absorb(Matrix& other_values, Vector& other_indices) 
    {
	const value_type    zero= math::zero(value_type());
#ifndef NDEBUG
	bool sorted = true;
	for(unsigned int i = 1; i < size(other_indices); ++i) {
	    sorted &= ( other_indices(i-1) < other_indices(i) );
	}
	assert(sorted);
#endif
		
	const std::size_t other_idx_size = size( other_indices );

	// Determine set of common indices.
	const int max_common_idx =
	    int(nb_vars() < other_idx_size ? nb_vars() : other_idx_size);
	mtl::dense_vector<int> my_idx( max_common_idx );
	mtl::dense_vector<int> ot_idx( max_common_idx );
	int offset = 0;
	for(int i = 0, j = 0; i < int(nb_vars()) && j < int(other_idx_size); ) {
	    int diff = (get_indices()(i) - other_indices(j));
	    if(diff == 0) {
		my_idx(offset) = i;
		ot_idx(offset) = j;
		++offset;
		++i;
		++j;
	    } else if(diff < 0) {
		++i;
	    } else {
		++j;
	    }
	}

	// Absorb the values.
	for(int i = 0; i < offset; ++i) {
	    for(int j = i; j < offset; ++j) {
		get_values()( my_idx(i), my_idx(j) ) += other_values( ot_idx(i), ot_idx(j) );
		other_values( ot_idx(i), ot_idx(j) ) = zero;
		get_values()( my_idx(j), my_idx(i) ) += other_values( ot_idx(j), ot_idx(i) );
		other_values( ot_idx(j), ot_idx(i) ) = zero;
	    }
	}
    }
    
    /// Removes the numerical values from the element.
    void clear() {
	m_neighbors.clear();
	m_neighbors.resize(1);
	matrix_type empty;
	swap(m_values, empty);
	m_indices.change_dim(0);
    }

    template< class ValueType > friend class element_structure;

  private:
    /// The set of neighbors of the element.
    neighbor_collection_type m_neighbors;

    /// The set of indices of this element.
    index_type m_indices;

    /// The [Size x Size] element matrix.
    matrix_type m_values;

    /// A unique sequence number for the element, indicating it's order relative to other elements.
    int m_sequence_number;

    int *dummy;
};


/// Print an element to an output stream.
template<typename OStream, class ValueType>
OStream& operator<<(OStream& out, element<ValueType>& el) 
{
    out << "ID: " << el.get_id() << "\n";
    if(el.nb_vars() > 0) {
	out << "Indices: (" << el.get_indices()(0);
	for(std::size_t i = 1; i < el.nb_vars(); ++i) {
	    out << ", " << el.get_indices()(i);
	}
	out << ")\n";
    } else {
	out << "Indices: ()\n";
    }
    out << "Neighbors: (";
    if (el.nb_vars() > 0) 
	for(int i = 0; i < el.get_nb_neighbors(); ++i) 
	    out << el.get_neighbors()[i]->get_id() << (i+1 < el.get_nb_neighbors()? ", " : ")\n");
    out << "Values: \n" << el.get_values();
    return out;
}

  }} // mtl::matrix

#endif // MTL_ELEMENT_INCLUDE
