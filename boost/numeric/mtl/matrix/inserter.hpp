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

#ifndef MTL_MATRIX_INSERTER_INCLUDE
#define MTL_MATRIX_INSERTER_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/operation/update.hpp>
#include <boost/numeric/mtl/detail/trivial_inserter.hpp>



namespace mtl { namespace mat {

/// Matrix inserter
/** The matrix inserter has two template arguments: the type of the target matrix and
    an update functor.
    The update functor determines how an existing entry is updated: overwritten, added,
    subtracted...
    The default is to overwrite existing entries.
**/
template <typename Matrix,
	  typename Updater = mtl::operations::update_store<typename Matrix::value_type> >
struct inserter 
  : public mtl::detail::trivial_inserter<Matrix, Updater>
{
    typedef mtl::detail::trivial_inserter<Matrix, Updater>     base;
    typedef typename Matrix::size_type                         size_type;

    explicit inserter(Matrix& matrix, size_type slot_size = 0) : base(matrix, slot_size) {}
};


template <typename Value, typename Parameters, typename Updater>
struct inserter<compressed2D<Value, Parameters>, Updater>
  : compressed2D_inserter<Value, Parameters, Updater>
{
    typedef compressed2D<Value, Parameters>                    matrix_type;
    typedef typename matrix_type::size_type                    size_type;
    typedef compressed2D_inserter<Value, Parameters, Updater > base;

    explicit inserter(matrix_type& matrix, size_type slot_size = 5) : base(matrix, slot_size) {}
};

template <typename Value, typename Parameters, typename Updater>
struct inserter<coordinate2D<Value, Parameters>, Updater>
  : coordinate2D_inserter<coordinate2D<Value, Parameters>, Updater>
{
    typedef coordinate2D<Value, Parameters>                                  matrix_type;
    typedef typename Parameters::size_type                                   size_type;
    typedef coordinate2D_inserter<coordinate2D<Value, Parameters>, Updater>  base;

    explicit inserter(matrix_type& matrix, size_type slot_size= 1) : base(matrix, slot_size) {}
};

template <typename Value, typename Parameters, typename Updater>
struct inserter<sparse_banded<Value, Parameters>, Updater>
  : sparse_banded_inserter<Value, Parameters, Updater>
{
    typedef sparse_banded<Value, Parameters>                                 matrix_type;
    typedef typename Parameters::size_type                                   size_type;
    typedef sparse_banded_inserter<Value, Parameters, Updater>  base;

    explicit inserter(matrix_type& matrix, size_type slot_size= 1) : base(matrix, slot_size) {}
};

template <typename Value, typename Parameters, typename Updater>
struct inserter<ell_matrix<Value, Parameters>, Updater>
  : ell_matrix_inserter<Value, Parameters, Updater>
{
    typedef ell_matrix<Value, Parameters>                    matrix_type;
    typedef typename matrix_type::size_type                    size_type;
    typedef ell_matrix_inserter<Value, Parameters, Updater > base;

    explicit inserter(matrix_type& matrix, size_type slot_size = 5) : base(matrix, slot_size) {}
};


}} // namespace mtl::matrix

#endif // MTL_MATRIX_INSERTER_INCLUDE
