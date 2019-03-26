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

#ifndef MTL_MATRICES_INCLUDE
#define MTL_MATRICES_INCLUDE

#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp> 
#include <boost/numeric/mtl/matrix/compressed2D.hpp> 
#include <boost/numeric/mtl/matrix/sparse_banded.hpp> 
#include <boost/numeric/mtl/matrix/multi_vector.hpp> 
#include <boost/numeric/mtl/matrix/multi_vector_range.hpp> 
#include <boost/numeric/mtl/matrix/element_matrix.hpp> 
#include <boost/numeric/mtl/matrix/element.hpp> 
#include <boost/numeric/mtl/matrix/implicit_dense.hpp> 
#include <boost/numeric/mtl/matrix/block_diagonal2D.hpp>
#include <boost/numeric/mtl/matrix/ell_matrix.hpp>
#include <boost/numeric/mtl/matrix/coordinate2D.hpp>

#include <boost/numeric/mtl/matrix/inserter.hpp> 
#include <boost/numeric/mtl/matrix/shifted_inserter.hpp> 
#include <boost/numeric/mtl/matrix/mapped_inserter.hpp>

#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/hermitian_view.hpp>
#include <boost/numeric/mtl/matrix/banded_view.hpp>
#include <boost/numeric/mtl/matrix/indirect.hpp>

#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/laplacian_setup.hpp> 
#include <boost/numeric/mtl/matrix/hessian_setup.hpp> 
#include <boost/numeric/mtl/matrix/poisson2D_dirichlet.hpp> 
#include <boost/numeric/mtl/matrix/identity2D.hpp> 

#include <boost/numeric/mtl/recursion/predefined_masks.hpp>

#include <boost/numeric/mtl/utility/matrix_type_generator.hpp>

#endif // MTL_MATRICES_INCLUDE
