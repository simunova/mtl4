
#ifndef ELEMENT_STRUCTURE_ALGORITHMS
#define ELEMENT_STRUCTURE_ALGORITHMS


#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/matrix/element.hpp>

namespace imf {

// Extracts an artificial element structure from the given sparse matrix. This
// algorithm
template< class Element_struct, class Matrix >
void greedy_extract_element_structure(Element_struct& es,  Matrix& M , std::string& output) {

/*******************************************************************************
 * PHASE 1: Generate Artificial 1x1 and kxk Elements
 ******************************************************************************/

	typedef unsigned int usint;

	// Elements
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type size_type;
	typedef mtl::mat::element<value_type> element_type;
	typedef typename element_type::index_type index_type;
	typedef typename element_type::matrix_type matrix_type;
	// typedef typename element_type::neighbor_iterator neigh_iterator;
				        
#if 0
	// Sparse matrices
	typedef glas::compressed_sparse_structure<glas::row_orientation>			sparse_structure;
	typedef typename sparse_structure::compressed_index_array_type			compressed_index_array_type;
	typedef glas::sparse_matrix<value_type, sparse_structure> sparse_type;

	compressed_index_array_type& row_start = glas::compressed_index_array( M );
	typename sparse_structure::index_array_type& col_idx =		glas::index_array( M );
	typename sparse_type::value_array_type& values = glas::value_array(M);
#endif
	std::vector<size_type> row_start(M.ref_major());
	std::vector<size_type> col_idx(M.ref_minor());
	std::vector<value_type> values(M.data);
	
	const value_type ZERO = value_type(0);
	const usint nb_rows = num_rows( M );

	// Determine a "suitable" guess for the maximum size of the elements.
	double avg_nnz_row = M.nnz();
	double std_dev = 0;
	for(usint r = 0; r < nb_rows; ++r) {
		std_dev += pow(double(row_start[r+1]-row_start[r]), 2) / double(nb_rows);
	}
	int max_idx = int(avg_nnz_row + 2*std_dev) + 1;
	if( max_idx < 7 ) {
		max_idx = 7;
	}

	std::cout << "max_idx="<<max_idx<< std::endl;

 	std::cout << "nnz="<< M.nnz() << std::endl;

	std::vector<element_type*> elements;
	elements.reserve(M.nnz() / max_idx);

	std::vector<int> idx;
	idx.reserve(max_idx);

	// Generate the element matrices from the given sparse matrix.
	int seq_nbr = 0;
	for(usint row = 0; row < nb_rows; ++row) {
		for(
			usint col_off = row_start[row];
				col_off < row_start[row+1];
		) {
			std::cout<< "row="<< row << "-"<< nb_rows << "\n";
			// We have at most k different indices.
			idx.clear();
			idx.push_back(row);
			for(
				int my_size = 1;
				(my_size < max_idx) && (col_off < row_start[row+1]);
				++col_off
			) {
				if( col_idx[col_off] != row ) {
					idx.push_back( col_idx[col_off] );
					++my_size;
				}
			}
			std::sort( idx.begin(), idx.end() );

			// Construct index.
			index_type el_index( idx.size() );
			for(usint i = 0; i < idx.size(); ++i) {
				el_index(i) = idx[i];
			}

			// Copy values.
			const usint n =  size(el_index);
			matrix_type el_vals(n,n);
			el_vals = ZERO;
			bool non_zero = false;
			for(usint i = 0; i < n; ++i) {
				const int lrow = el_index(i);

				// Seek matching indices.
				usint k = row_start[lrow];
				usint j = 0;
				while( (k < row_start[lrow+1]) && (j < n) ) {
					const int diff = el_index(j) - col_idx[k];
					if(diff > 0) {
						++k;
					} else if(diff < 0) {
						++j;
					} else {
						// We found a match. Copy the value.
						el_vals( i,j ) = values[k];
						values[k] = ZERO;
						non_zero |= (el_vals(i,j) != ZERO);
						++j;
						++k;
					}
				}
			}

			// if the element matrix is non-zero, add it.
			if(non_zero) {
				element_type* el = new element_type(seq_nbr, el_index, el_vals);
				elements.push_back(el);
				++seq_nbr;
			}
		}
		// If the diagonal entry is not included in the sparsity of the matrix,
		// it is wise *not* to add an element with a single node. In this
		// manner, we will avoid the situation wherein that element is selected,
		// and the matrix is zero at this position, but where the elimination of
		// a different element (which includes the node of this row as off-
		// diagonal element) would make the matrix at this position non-zero
		// (and hence invertible).
	}

	// Remove unnecessary nodes.
	for(usint i = 0; i < elements.size(); ++i) {
		element_type& el = *(elements[i]);
		index_type& idx = el.get_indices();
		matrix_type& vals = el.get_values();

		std::vector<int> remove_nodes;
		for(int c = 0; c < el.nb_vars(); ++c) {
			bool all_zero = true;
			for(int k = 0; k < el.nb_vars(); ++k) {
				all_zero &= ( (vals(k,c) == ZERO) && (vals(c,k) == ZERO) );
			}
			if(all_zero) {
				remove_nodes.push_back( idx(c) );
			}
		}
		el.remove_nodes(remove_nodes, el);
	}

	// Generate neighbourhood information.
	std::vector<int>* node_element_map = new std::vector<int>[nb_rows];
	for( usint i = 0; i < elements.size(); ++i ) {
		element_type& el = *(elements[i]);
		index_type& idx = el.get_indices();
		for(int j = 0; j < el.nb_vars(); ++j) {
			node_element_map[ idx(j) ].push_back(el.get_id());
		}
	}
	// Construct neighbourhood set for every element.
	for( usint i = 0; i < elements.size(); ++i ) {
		element_type& el = *(elements[i]);
		index_type& idx = el.get_indices();
		std::set<int> neighs;
		for(int j = 0; j < el.nb_vars(); ++j) {
			neighs.insert(
				node_element_map[ idx(j) ].begin(),
				node_element_map[ idx(j) ].end()
			);
		}
			// An element shouldn't neighbour itself.
		neighs.erase( el.get_id() );
		for(
			std::set<int>::iterator it = neighs.begin();
			it != neighs.end();
			++it
		) {
			el.get_neighbors().push_back( elements[*it] );
		}
	}
	delete[] node_element_map;

/*******************************************************************************
 * PHASE 3: Construct Grid
 ******************************************************************************/
	const int orig_els = elements.size();
	element_type* elems = new element_type[ elements.size() ];

	// Map the location of the new elements.
	std::vector<int> seq_nbr_map( orig_els, -1 );
	for(usint i = 0; i < elements.size(); ++i) {
		seq_nbr_map[ elements[i]->get_id() ] = i;
	}

	// Construct the element set.
	for(usint i = 0; i < elements.size(); ++i) {
		// Copy the element to the grid.
		elems[i] = *(elements[i]);
		element_type& el = elems[i];
		el.get_id() = i;

		typename element_type::neighbor_collection_type& neighs =
			el.get_neighbors();
		for( usint j = 0; j < neighs.size(); ++j ) {
			// Remap the neighbour pointer to the new location.
			int offset = seq_nbr_map[neighs[j]->get_id()];
			neighs[j] = elems + offset;

			assert( offset >= 0 );
		}

		// Clear the data of the copy.
		elements[i]->clear();
	}
#if 0
#ifndef NDEBUG
#if VERBOSE_MODE > 10
	std::cout << "The elements in the algebraic grid: " << std::endl;
	for(usint i = 0; i < elements.size(); ++i) {
		std::cout << elems[i] << std::endl;
	}
#endif
#endif
#endif

	es.consume(orig_els, nb_rows, elems);
#if 0
	element_structure<value_type>* es = new element_structure<value_type>(
		orig_els, nb_rows, elems
	);
#endif
	// Delete the copied elements.
	for(usint i = 0; i < elements.size(); ++i) {
		assert( elements[i] );
		delete elements[i];
	}
	elements.clear();

	std::cout << "Generated " << orig_els << " artificial elements." <<
		std::endl;
 es.write_to_file(output);
//	return es;
}

} // end namespace imf

#endif /* ELEMENT_STRUCTURE_ALGORITHMS */
