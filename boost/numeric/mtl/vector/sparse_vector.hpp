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

#ifndef MTL_VECTOR_SPARSE_VECTOR_INCLUDE
#define MTL_VECTOR_SPARSE_VECTOR_INCLUDE

#include <iostream>
#include <vector>
#include <algorithm>

#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/zipped_sort.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace vec {

/// This is not (yet) a full numeric class!!! It is only a helper for factorization with very limited functionality. DO NOT USE in numeric code (yet)!!!
template <typename Value, typename Parameter= parameters<> >
class sparse_vector
{
    typedef Value                           value_type;
    typedef std::vector<value_type>         vv_type;
    typedef typename Parameter::size_type   size_type;
    typedef std::vector<size_type>          sv_type;
    typedef typename Parameter::orientation orientation;
    typedef sparse_vector                   self;

    void check(size_type MTL_DEBUG_ARG(j)) const
    {   MTL_DEBUG_THROW_IF(is_negative(j) || j > size_type(indices.size()), mtl::index_out_of_range()); }

  public:
    sparse_vector(size_type n= 0) : n(n), on_indices(true) /*, longest(0) */ {}

    // ~sparse_vector() { std::cout << "longest vector was " << longest << '\n'; }

    std::size_t nnz() const { assert(indices.size() == data.size()); return indices.size(); }

    bool exists(size_type i) const
    {   return find(indices.begin(), indices.end(), i) != indices.end(); }

    value_type operator[](size_type i) const
    {
	typename sv_type::iterator it= find(indices.begin(), indices.end(), i);
	// size_type j= distance(indices.begin(), it);
	return it != indices.end() ? data[distance(indices.begin(), it)] : value_type(0);
    }

    void insert(size_type i, const value_type& v)
    {
	// vampir_trace<9901> tracer;
	typename sv_type::iterator it= std::lower_bound(indices.begin(), indices.end(), i);
	if (it == indices.end()) {
	    // std::cout << "Insert entry " << i << " at the end\n";
	    indices.push_back(i);
	    data.push_back(v);
	    return;
	}

	// std::cout << "Insert entry " << i << " at position " << distance(indices.begin(), it) << '\n';
	typename vv_type::iterator dit= data.begin(); 
	// std::cout << "*dit == " << *dit << '\n';
	advance(dit, distance(indices.begin(), it));
	// std::cout << "*dit == " << *dit << ", distance(data.begin(), dit) == " << distance(data.begin(), dit) << '\n';
	data.insert(dit, v); // insert v at same position
	indices.insert(it, i);
    }

    size_type pos(size_type i) const 
    {
	typename sv_type::const_iterator it= std::lower_bound(indices.begin(), indices.end(), i);
	MTL_DEBUG_THROW_IF(it == indices.end(), logic_error("Position doesn't exists"));
	return distance(indices.begin(), it);
    }

    value_type& operator[](size_type i)
    {
	if (!exists(i)) 
	    insert(i, value_type(0));
	return data[pos(i)];
    }

    void make_empty() { data.resize(0); indices.resize(0); }

    void crop(value_type threshold= value_type(0))
    {
	using std::abs;
	std::size_t j= 0;
	for (std::size_t i= 0, end= data.size(); i < end; i++) 
	    if (abs(data[i]) > threshold) {
		if (i > j) {
		    data[j]= data[i];
		    indices[j]= indices[i];
		}
		++j;
	    }
	data.resize(j);
	indices.resize(j);
    }

    void sort_on_data() // largest magnitude first
    {
	// vampir_trace<9903> tracer;
	if (n > 0) {

	    std::sort(utility::zip_it<value_type, size_type>(&data[0], &indices[0], 0),
		utility::zip_it<value_type, size_type>(&data[0], &indices[0], indices.size()),
		utility::abs_greater_0());
	}
	on_indices= false;

	// if (indices.size() > longest) longest= indices.size();
    }

    void sort_on_indices()
    {
	if (n > 0)
	    std::sort(utility::zip_it<size_type, value_type>(&indices[0], &data[0], 0), 
		      utility::zip_it<size_type, value_type>(&indices[0], &data[0], indices.size()), 
		      utility::less_0());
	on_indices= true;
    }

    size_type index(size_type j) const { check(j); return indices[j]; }
    value_type const& value(size_type j) const { check(j); return data[j]; }
    value_type&       value(size_type j)       { check(j); return data[j]; }
    
    std::pair<size_type, value_type> entry(size_type j) const
    {   return std::make_pair(indices[j], data[j]);    }

    std::pair<size_type, value_type> largest(size_type j) const
    {
	MTL_DEBUG_THROW_IF(on_indices, logic_error("Vector must be sorted on values to use this function"));
	MTL_DEBUG_THROW_IF(j > size_type(data.size()), logic_error("There are less than j entries."));
	return std::make_pair(indices[j], data[j]);
    }

    friend std::ostream& operator<<(std::ostream& out, const self& v)
    {
	out << '{' << v.n << (traits::is_row_major<orientation>::value ? "R" : "C") << "}[";
	for (std::size_t i= 0, end= v.data.size(); i < end; i++) 
	    out << v.indices[i] << ':' << v.data[i] << (i+1 < end ? ", " : "");
	return out << ']';
    }

  private:
    size_type   n;
    sv_type     indices;
    vv_type     data;
    bool        on_indices; // if true sorted on indices, otherwise on values
    
    // std::size_t longest;
};

template <typename Value, typename Parameter>
std::size_t inline size(const sparse_vector<Value, Parameter>& v)
{   return v.n; }


}} // namespace mtl::vector

namespace mtl {
    using vec::sparse_vector;
}

#endif // MTL_VECTOR_SPARSE_VECTOR_INCLUDE
