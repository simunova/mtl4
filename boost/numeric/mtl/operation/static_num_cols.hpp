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

#ifndef MTL_STATIC_NUM_COLS_INCLUDE
#define MTL_STATIC_NUM_COLS_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/matrix/dimension.hpp>
#include <boost/numeric/mtl/vector/dimension.hpp>

namespace mtl {

/// Number of columns given at compile time
/** General declaration, used to disable unsupported types **/
template <typename Collection>
struct static_num_cols {
    // typedef xxx type;
    // static const type value= yyy;
};
    

/// static_num_cols implementation for (1D) arrays interpreted as vectors
template <typename Value, unsigned Size>
struct static_num_cols<Value[Size]>
{
    typedef std::size_t   type;
    static const type value= 1;
};	   

/// static_num_cols implementation for (2D and higher) arrays interpreted as matrices
template <typename Value, unsigned Rows, unsigned Cols>
struct static_num_cols<Value[Rows][Cols]>
{
    typedef std::size_t   type;
    static const type value= Cols;    
};	    


template <std::size_t Size>
struct static_num_cols< vec::fixed::dimension<Size> >
{
    typedef std::size_t   type;
    static const type value= Size;    
};

template <typename V, typename P> 
struct static_num_cols<mtl::vec::dense_vector<V, P> > 
{
    typedef std::size_t   type;
    static const type value= traits::is_row_major<P>::value ? static_num_cols<typename P::dimension>::value : 1;   
};


template <std::size_t Rows, std::size_t Cols>
struct static_num_cols< fixed::dimensions<Rows, Cols> >
{
    typedef std::size_t   type;
    static const type value= Cols;    
};

template <typename V, typename P> 
struct static_num_cols<mtl::mat::dense2D<V, P> > 
  : static_num_cols<typename P::dimensions> 
{};

template <typename V, std::size_t M, typename P> 
struct static_num_cols<mtl::mat::morton_dense<V, M, P> > 
  : static_num_cols<typename P::dimensions> 
{};

template <typename V, typename P> 
struct static_num_cols<mtl::mat::compressed2D<V, P> > 
  : static_num_cols<typename P::dimensions> 
{};

} // namespace mtl

#endif // MTL_STATIC_NUM_COLS_INCLUDE
