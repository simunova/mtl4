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



#ifndef MTL_ELEMENT_STRUCTURE_INCLUDE
#define MTL_ELEMENT_STRUCTURE_INCLUDE

#include <iostream>
#include <ostream>

#include <boost/numeric/mtl/matrix/element.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>

namespace mtl { namespace mat {

#if 0   ///TODO need for write elmement_structur
namespace print {

	template< class Type >
	struct print_type {
		template<class Stream>
		static void print(Stream&);
	};

	template< >
	struct print_type<double> {
		template<class Stream>
		static void print(Stream& str) {
			str << "double\n";
		}
	};

	template< >
	struct print_type<float> {
		template<class Stream>
		static void print(Stream& str) {
			str << "double\n";
		}
	};

	template< >
	struct print_type<std::complex<double> > {
		template<class Stream>
		static void print(Stream& str) {
			str << "complex\n";
		}
	};

	template< class Type >
	struct print_value {
		template<class Stream>
		static void print(Stream&);
	};

	template< >
	struct print_value<double> {
		template<class Stream>
		static void print(Stream& str, double& val) {
			str << std::scientific << std::setprecision(15);
			str << val;
		}
	};

	template< >
	struct print_value<float> {
		template<class Stream>
		static void print(Stream& str, float& val) {
			str << std::scientific << std::setprecision(15);
			str << double(val);
		}
	};

	template< >
	struct print_value<std::complex<double> > {
		template<class Stream>
		static void print(Stream& str, std::complex<double>& val) {
			str << std::scientific << std::setprecision(15);
			str << val.real() << "\t" << val.imag() << "\t";
		}
	};
}
  #endif

/**
 * A generic abstract base class for meshes. It describes the concept of a mesh.
 */
template< class ValueType >
class element_structure 
{

public:
    /// Type of the numerical values of the element coefficient matrices.
    typedef ValueType value_type;

    /// Type of the element.
    typedef element<value_type> element_type;

    /// Type of index arrays.
    typedef typename element_type::index_type index_type;
    // Type of the indices themselves
    typedef typename Collection<index_type>::value_type ii_type;

    /// Type of the iterator over the elements of the mesh.
    typedef element_type* element_iterator;

    /// Type of this class.
    typedef element_structure<ValueType> this_type;
    typedef this_type                    self;
    
    /// Standard constructor.
    element_structure(int total_elements= 0, int total_vars= 0, element_type* elements= 0)
      : m_total_elements(total_elements), m_total_vars(total_vars),
	m_elements(elements), index_heap(0), value_heap(0)
    { }

    /// consume elements into element_structure
    void consume(int total_elements, int total_vars, element_type* elements)
    {
	m_total_elements= total_elements;
	m_total_vars= total_vars;
	delete[] m_elements;
	m_elements= elements;
	delete[] index_heap; index_heap= 0;
	delete[] value_heap; value_heap= 0;
    }


    /// Copy the given mesh.
    element_structure(this_type const& other)
      : m_total_elements(other.m_total_elements),
	m_total_vars(other.m_total_vars),
	m_elements(m_total_elements == 0 ? 0 : new element_type[m_total_elements]), 
	index_heap(0), value_heap(0)
    {
	typedef typename element_type::neighbor_collection_type neigh_coll_type;

	int j = 0;
	bool ordered = true;
	for(element_iterator it = other.element_begin(); it != other.element_end(); ++it) {
	    // Deep copy the elements.
	    m_elements[j] = *it;
	    ordered &= (it->get_id() == j);
	    ++j;
	}
	assert( ordered );
	// Reconstruct the network of neighbors.
	for(element_iterator it = this->element_begin(); it != this->element_end(); ++it) {
	    neigh_coll_type new_neighs;
	    neigh_coll_type& old_neighs = it->get_neighbors();
	    for(int i = 0; i < it->get_nb_neighbors(); ++i) {
		element_type& neigh = *(old_neighs[i]);
		int pos = neigh.get_id();
		new_neighs.push_back( this->m_elements+pos );
	    }
	    old_neighs.assign(new_neighs.begin(), new_neighs.end());
	}
    }

    /// Destructor
    ~element_structure() { delete[] m_elements; delete[] index_heap; delete[] value_heap; }

    /// make compakt memory block from elements
    void make_compact()
    {
	assert(index_heap == 0); assert(value_heap == 0); // might be relaxed later
	
	std::size_t total_indices= 0, total_values= 0;
	for (int i= 0; i < m_total_elements; i++) {
	    total_indices+= m_elements[i].nb_vars();
	    total_values+= m_elements[i].nb_values();
	}
	index_heap= new ii_type[total_indices];
	value_heap= new value_type[total_values];

	int index_pos= 0, value_pos= 0;
	for (int i= 0; i < m_total_elements; i++) {
	    element_type& element= m_elements[i];
	    int s= int(element.nb_vars());
	    index_type index_tmp(s, index_heap + index_pos);
	    index_tmp= element.get_indices();
	    swap(index_tmp, element.get_indices());
	    index_pos+= s;

	    typename element_type::matrix_type value_tmp(s, s, value_heap + value_pos);
	    value_tmp= element.get_values();
	    swap(value_tmp, element.get_values());
	    value_pos+= s * s;
	}
	assert(total_indices == std::size_t(index_pos));
	assert(total_values  == std::size_t(value_pos));
    }

    /*******************************************************************************
     * Inspector Members
     ******************************************************************************/

    /// Total number of elements in the grid.
    int get_total_elements() const { return m_total_elements; }

    /// Total number of variables.
    int get_total_vars() const { return m_total_vars;   }

    /// Total number of non-zero values.
    int get_total_nnz() const 
    {
	int nnz = 0;
	for(element_iterator it = element_begin(); it != element_end(); ++it) {
	    nnz += it->nnz();
	}
	return nnz;
    }

    /// Iterator to the first element.
    element_iterator element_begin() const { return m_elements + 0;  }

    /// An iterator to the element past the last element.
    element_iterator element_end() const { return m_elements + this->get_total_elements();   }

#if 1
    /// Writes the elements to the specified file.  TODO at the moment very slow
    void write_to_file(const std::string& filename) 
    {
	//using namespace print;

	std::ofstream file(filename.c_str());

	// Write header information.
	file << get_total_elements() << "\n";
	file << this->get_total_vars() << "\n";
	//print_type<value_type>::print(file);

	// Write element matrices.
	for(element_iterator it = element_begin(); it != element_end(); ++it) {
	    // Write indices.
	    for(int i = 0; i < it->nb_vars()-1; ++i) {
		file << it->get_indices()(i) << " ";
	    }
	    file << it->get_indices()(it->nb_vars()-1) << "\n";

	    // Write values.
	    for(int r = 0; r < it->nb_vars(); ++r) {
		for(int c = 0; c < it->nb_vars()-1; ++c) {
		    //print_value<value_type>::print(file, it->get_values()(r,c));
		    file << it->get_values()(r,c); file << " ";
		}
		//print_value<value_type>::print(file, it->get_values()(r,it->nb_vars()-1));
		file << it->get_values()(r,it->nb_vars()-1);
		file << "\n";
	    }
	    file << "\n";
	}
    }
#endif
        
    int           m_total_elements; ///< The total number of elements.
    int           m_total_vars; ///< The total number of variables.
    element_type* m_elements; ///< The elements of the grid, stored consecutively.

    ii_type*      index_heap;
    value_type*   value_heap;
};

template <typename ValueType>
inline std::size_t num_rows(const element_structure<ValueType>& A)
{   return A.get_total_vars(); }

template <typename ValueType>
inline std::size_t num_cols(const element_structure<ValueType>& A)
{   return A.get_total_vars(); }

template <typename ValueType>
inline std::size_t size(const element_structure<ValueType>& A)
{   return A.get_total_vars() * A.get_total_vars(); }


template <typename ValueType>
inline void swap(element_structure<ValueType>& x, element_structure<ValueType>& y)
{
    swap(x.m_total_elements, y.m_total_elements);
    swap(x.m_total_vars, y.m_total_vars);
    swap(x.m_elements, y.m_elements);
}


}} // mtl::matrix


#endif // MTL_ELEMENT_STRUCTURE_INCLUDE
