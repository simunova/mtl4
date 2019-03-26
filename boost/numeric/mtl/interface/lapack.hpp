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

#ifndef MTL_LAPACK_INCLUDE
#define MTL_LAPACK_INCLUDE

#include <complex>


#ifdef __cplusplus
extern "C" {
#endif


  // Cholesky Factorization
  void spotrf_(const char* uplo, const int* n, float *a, const int* ld, int* info);
  void dpotrf_(const char* uplo, const int* n, double *a, const int* ld, int* info);
  void cpotrf_(const char* uplo, const int* n, std::complex<float> *a, const int* ld, int* info);
  void zpotrf_(const char* uplo, const int* n, std::complex<double> *a, const int* ld, int* info);














#ifdef __cplusplus
} // extern "C"
#endif


#endif // MTL_LAPACK_INCLUDE
