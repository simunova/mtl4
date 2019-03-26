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

#ifndef ITL_REPEATED_INCLUDE
#define ITL_REPEATED_INCLUDE

namespace itl {

/// Repeat the smoother
template <typename Smoother, std::size_t N= 1>
class repeated
{
  public:

    typedef Smoother smoother_type; 

    /// Construct with \p smoother
    repeated(const Smoother& smoother) : smoother(smoother) {}

    /// Construct with \p smoother and number of repetitions \p n
    template <typename Matrix>
    repeated(const Matrix& A) : smoother(A) {}

    /// Apply smoother n times on vector \p x, i.e. \p x is changed
    template <typename Vector, typename RHSVector>
    Vector& operator()(Vector& x, const RHSVector& b) const
    {
	for (std::size_t i= 0; i < N; i++)
	    smoother(x, b);
 	return x;
    }

  private:
    Smoother    smoother;
    std::size_t n;
};


} // namespace itl

#endif // ITL_REPEATED_INCLUDE
