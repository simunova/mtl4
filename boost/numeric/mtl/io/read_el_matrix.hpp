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

#ifndef MTL_MATRIX_READ_EL_MATRIX
#define MTL_MATRIX_READ_EL_MATRIX

#include <string>

#include <iostream>
#include <istream>
#include <set>
#include <vector>
#include <valarray>

#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/mtl/matrix/element.hpp>
#include <boost/numeric/mtl/matrix/element_structure.hpp>

namespace mtl { namespace mat {

// Read a value from the stream. The stream is advanced.
template <class T, class StreamType>
inline T read_value(StreamType& stream) 
{
	T value;
	stream >> value;
	return value;
}

// Reads the element structure from a given file.
//
// It is assumed the nodes are numbered consecutively, i.e. there are no unused
// node numbers.
template < typename StreamType, typename ValueType>
void read_el_matrix(StreamType& file, element_structure<ValueType>& A) 
{
    // Type definitions
    typedef element<ValueType>		element_type;
    // typedef typename element_type::value_type 	value_type;
    typedef typename element_type::index_type 	indices;
    typedef typename element_type::matrix_type 	matrix;
    vampir_trace<4036> trace;

    // Read element type information.
    int nb_elements = 0;
    file >> nb_elements;
    file.ignore(500,'\n');
    std::cout << "nb elements: " << nb_elements << "\n";

    // Compatibility with older files
    file.ignore(500,'\n');

    assert(nb_elements >= 0);

    element_type* elements = new element_type[nb_elements];

    // Read elements from file.
    int el_nbr = 0;
    int nb_total_vars = 0;
    while( el_nbr < nb_elements ) {

	// Read the node numbers.
	std::string line;
	getline(file, line, '\n');
	std::stringstream node_line(line), read_node_line(line);

	int read_num=0, i=0;
	while( !read_node_line.eof() ) {
	    int idx = 0;
	    read_node_line >> idx;
	    ++read_num;
	}
	read_num--;
	mtl::dense_vector<int> nodes(read_num, 0);  
	while( !node_line.eof() ) {
	    int idx = 0;
	    node_line >> idx;
	    if (i<read_num)
		nodes[i]=idx;
	    if(idx > nb_total_vars) 
		nb_total_vars = idx;	    
	    i++;
	}
	indices index(nodes);
	// Read the values.
	const int nb_vars = int(size(nodes));
	matrix vals(nb_vars, nb_vars);
	for(int i = 0; i < nb_vars*nb_vars; ++i) 
	    vals(i / nb_vars, i % nb_vars) = read_value<ValueType>(file);	
	file.ignore(500,'\n');
	file.ignore(500,'\n');
	element_type elem(el_nbr, index, vals);
	elements[el_nbr] = elem;
	if(el_nbr == 0){
	  std::cout<< "elem=" << elem << "\n";
	}
	++el_nbr;
    }

    // Construct mapping.
    ++nb_total_vars;
    assert(nb_total_vars >= 0);
    std::vector<int>* node_element_map = new std::vector<int>[nb_total_vars];
    for( int i = 0; i < nb_elements; ++i ) {
	element_type& el = elements[i];
	indices& idx = el.get_indices();
	for(std::size_t j = 0; j < el.nb_vars(); ++j) 
	    node_element_map[ idx(j) ].push_back(el.get_id());	
    }

    // Construct neighborhood information.
    for( int i = 0; i < nb_elements; ++i ) {
	element_type& el = elements[i];
	indices& idx = el.get_indices();
	std::set<int> neighs;
	for(std::size_t j = 0; j < el.nb_vars(); ++j) 
	    neighs.insert(node_element_map[ idx(j) ].begin(),
			  node_element_map[ idx(j) ].end());
	
	for(std::set<int>::iterator it = neighs.begin(); it != neighs.end(); ++it) 
	    if( *it != el.get_id() ) 
		el.get_neighbors().push_back( elements+(*it) );
	    	

	// Sort data.
	el.sort_indices();
    }

    delete[] node_element_map;
    A.consume(nb_elements, nb_total_vars, elements);
}

template <typename ValueType>
inline void read_el_matrix(std::string& mat_file, element_structure<ValueType>& A) 
{    read_el_matrix(mat_file.c_str(), A);   }

template <typename ValueType>
void read_el_matrix(const char* mat_file, element_structure<ValueType>& A) 
{
    std::ifstream file;
    file.open( mat_file );
    if( !file.is_open() ) {
	std::cout << "The file \"" << mat_file << "\" could not be opened." <<
	    std::endl;
	throw "File could not be opened"; 
    }
    read_el_matrix(file, A);

    file.close();
}

}} // end namespace mtl::matrix

#endif // MTL_MATRIX_READ_EL_MATRIX
