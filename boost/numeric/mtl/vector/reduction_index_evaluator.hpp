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

#ifndef MTL_VECTOR_REDUCTION_INDEX_EVALUATOR_INCLUDE
#define MTL_VECTOR_REDUCTION_INDEX_EVALUATOR_INCLUDE

namespace mtl { namespace vec {

/// Class for index-wise vector reductions
template <typename Scalar, typename Vector, typename Functor, typename Assign>
struct reduction_index_evaluator
{
    reduction_index_evaluator(Scalar& scalar, const Vector& v) 
      : scalar(scalar), v(v) 
    {
	Functor::init(tmp[0]);
	tmp[1]= tmp[2]= tmp[3]= tmp[0];
    }

    ~reduction_index_evaluator() 
    { 
	Functor::finish(tmp[0], tmp[1]);
	Functor::finish(tmp[2], tmp[3]);
	Functor::finish(tmp[0], tmp[2]);
	Assign::apply(scalar, Functor::post_reduction(tmp[0])); // compute sqrt or such if necessary
    }

    template <unsigned Offset>
    void at(std::size_t i) 
    { 
	Functor::update(tmp[Offset], v[i+Offset]); 
    }

    void operator[] (std::size_t i) { at<0>(i); }
    void operator() (std::size_t i) { at<0>(i); }    

    Scalar&        scalar;
    Scalar         tmp[4];
    const Vector&  v;
};

template <typename Scalar, typename Vector, typename Functor, typename Assign>
inline std::size_t size(const reduction_index_evaluator<Scalar, Vector, Functor, Assign>& eval)
{ return size(eval.v); }


}} // namespace mtl::vector

#endif // MTL_VECTOR_REDUCTION_INDEX_EVALUATOR_INCLUDE
