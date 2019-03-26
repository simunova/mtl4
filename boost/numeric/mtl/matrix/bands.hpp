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

#ifndef MTL_MATRIX_BANDS_INCLUDE
#define MTL_MATRIX_BANDS_INCLUDE

#include <boost/numeric/mtl/matrix/banded_view.hpp>

namespace mtl { namespace mat {

namespace traits {

    template <typename Matrix>
    struct bands
    {
	typedef banded_view<Matrix> type;
    };
}

/// Returns a view of a matrix \p A from diagonal \p begin to \p end
/** The main diagonal is numbered 0; the off-diagonal below the main one is -1.
    Accordingly, the off-diagonal above the main is 1.
    The parameters \p begin and \p end specify a right-open interval.
    For, instance bands(A, -1, 2) yields a tridiagonal matrix. **/
template <typename Matrix> 
typename traits::bands<Matrix>::type
inline bands(const Matrix& A, long begin, long end)
{
    typedef typename traits::bands<Matrix>::type result;
    return result(A, begin, end);
}

} // namespace matrix

    using mat::bands;

} // namespace mtl

#endif // MTL_MATRIX_BANDS_INCLUDE
