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

#ifndef MTL_VECTOR_DOT_INDEX_EVALUATOR_INCLUDE
#define MTL_VECTOR_DOT_INDEX_EVALUATOR_INCLUDE

namespace mtl { namespace vec {

/// Class for index-wise computation of dot product
template <typename Scalar, typename Vector1, typename Vector2, typename ConjOpt, typename Assign>
struct dot_index_evaluator
{
    dot_index_evaluator(Scalar& scalar, const Vector1& v1, const Vector2& v2) 
      : scalar(scalar), v1(v1), v2(v2) 
    { 
	tmp[0]= tmp[1]= tmp[2]= tmp[3]= Scalar(0); 
    }

    ~dot_index_evaluator() 
    { 
	Scalar s(tmp[0] + tmp[1] + tmp[2] + tmp[3]);
	Assign::apply(scalar, s); 
    }
    
    void operator() (std::size_t i) { tmp[0]+= ConjOpt()(v1[i]) * v2[i]; }
    void operator[] (std::size_t i) { (*this)(i); }

    template <unsigned Offset>
    void at(std::size_t i) 
    { tmp[Offset]+= ConjOpt()(v1[i+Offset]) * v2[i+Offset]; }

    Scalar&        scalar;
    Scalar         tmp[4];
    const Vector1& v1;
    const Vector2& v2;
};

template <typename Scalar, typename Vector1, typename Vector2, typename ConjOpt, typename Assign>
inline std::size_t size(const dot_index_evaluator<Scalar, Vector1, Vector2, ConjOpt, Assign>& eval)
{ 
    return size(eval.v1);
}


}} // namespace mtl::vector

#endif // MTL_VECTOR_DOT_INDEX_EVALUATOR_INCLUDE
