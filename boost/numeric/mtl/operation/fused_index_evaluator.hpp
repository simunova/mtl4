// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_FUSED_INDEX_EVALUATOR_INCLUDE
#define MTL_FUSED_INDEX_EVALUATOR_INCLUDE

#include <boost/numeric/mtl/utility/index_evaluator.hpp>

namespace mtl { namespace vec {

template <typename T, typename U>
struct fused_index_evaluator
{
    fused_index_evaluator(T& first, U& second) 
      : first(index_evaluator(first)), second(index_evaluator(second)) {}

    template <unsigned Offset>
    void at(std::size_t i) 
    {
	first.template at<Offset>(i);
	second.template at<Offset>(i);
    }

    void operator() (std::size_t i) { at<0>(i); }
    void operator[] (std::size_t i) { at<0>(i); }

    typename traits::index_evaluator<T>::type first;
    typename traits::index_evaluator<U>::type second;
};

template <typename T, typename U>
inline size_t size(const fused_index_evaluator<T, U>& expr) { return size(expr.first); }

}} // namespace mtl::vector

#endif // MTL_FUSED_INDEX_EVALUATOR_INCLUDE
