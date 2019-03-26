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


#ifndef MTL_BLAS_INCLUDE
#define MTL_BLAS_INCLUDE


#ifdef MTL_HAS_BLAS

// #include "mtl/mtl_config.h"

#if 0
#include "mtl/mtl_complex.h"
using std::complex;
#endif

/*--------------------------------------------------------
   Basic Linear Algebra Subprograms for C/C++
   Version 1.0
   Matthew E. Gaston
   May 6, 1998
----------------------------------------------------------*/

#define MTL_BLAS_NAME(name) name##_

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------
    Level 1 BLAS
-----------------------------------------------------------*/

/*
//  Dot product functions
*/
float MTL_BLAS_NAME(sdot)(int*, float*, int*, float*, int*);
double MTL_BLAS_NAME(dsdot)(int*, float*, int*, float*, int*);
float MTL_BLAS_NAME(sdsdot)(int*, float*, float*, int*, float*, int*);
double MTL_BLAS_NAME(ddot)(int*, double*, int*, double*, int*);

/*
    AXPY
*/
void MTL_BLAS_NAME(saxpy)(int*, float*, float*, int*, float*, int*);
void MTL_BLAS_NAME(daxpy)(int*, double*, double*, int*, double*, int*);

/*
    Copy
*/
void MTL_BLAS_NAME(scopy)(int*, float*, int*, float*, int*);
void MTL_BLAS_NAME(dcopy)(int*, double*, int*, double*, int*);

/*
    Swap
*/
void MTL_BLAS_NAME(sswap)(int*, float*, int*, float*, int*);
void MTL_BLAS_NAME(dswap)(int*, double*, int*, double*, int*);

/*
    2 Norm
*/
float MTL_BLAS_NAME(snrm2)(int *, float*, int*);
double MTL_BLAS_NAME(dnrm2)(int *, double*, int*);

/*
    Sum of Absolute Values
*/
float MTL_BLAS_NAME(sasum)(int *, float*, int*);
double MTL_BLAS_NAME(dasum)(int *, double*, int*);

/*
    Scale
*/
void MTL_BLAS_NAME(sscal)(int*, float*, float*, int*);
void MTL_BLAS_NAME(dscal)(int*, double*, double*, int*);

/*
    Maximum absolute value
*/
int MTL_BLAS_NAME(isamax)(int *, float*, int*);
int MTL_BLAS_NAME(idamax)(int *, double*, int*);


/*
    Givens Plane Rotation
*/
void MTL_BLAS_NAME(srotg)(float*, float*, float*, float*);
void MTL_BLAS_NAME(drotg)(double*, double*, double*, double*);
#if 0
void MTL_BLAS_NAME(crotg)(complex<float>*,complex<float>*,float*,complex<float>*);
void MTL_BLAS_NAME(zrotg)(complex<double>*,complex<double>*,double*,complex<double>*);
#endif
void MTL_BLAS_NAME(srot)(int*, float*, int*, float*, int*, float*, float*);
void MTL_BLAS_NAME(drot)(int*, double*, int*, double*, int*, double*, double*);
#if 0
/* MTL implements ccrot and zzrot */
void MTL_BLAS_NAME(csrot)(int*, complex<float>*, int*, complex<float>*, int*, 
	    complex<float>*, complex<float>*);
void MTL_BLAS_NAME(zdrot)(int*, complex<double>*, int*, complex<double>*, int*,
	    double*, double*);
#endif

/*---------------------------------------------------------
    Level 2 BLAS
-----------------------------------------------------------*/

void MTL_BLAS_NAME(dgemv)(char*, int*, int*, double*, double*, int*,
	    double*, int*, double*, double*, int*);

void MTL_BLAS_NAME(dger)(int*, int*, double*, double*, int*, double*,
	   int*, double*, int*);

void MTL_BLAS_NAME(dgbmv)(char*, int*, int*, int*, int*, double*, double*, int*,
	    double*, int*, double*, double*, int*);


void MTL_BLAS_NAME(dtrsv)(char* uplo, char* trans, char* diag, int* n, double *da, 
	    int* lda, double *dx, int* incx);

/*---------------------------------------------------------
    Level 3 BLAS
-----------------------------------------------------------*/

void MTL_BLAS_NAME(sgemm)(const char* transa, const char* transb, 
	    const int* m, const int* n, const int* k,
	    const float* alpha,  const float *da,  const int* lda,
	    const float *db, const int* ldb, const float* dbeta,
	    float *dc, const int* ldc);

void MTL_BLAS_NAME(dgemm)(const char* transa, const char* transb, 
	    const int* m, const int* n, const int* k,
	    const double* alpha,  const double *da,  const int* lda,
	    const double *db, const int* ldb, const double* dbeta,
	    double *dc, const int* ldc);

void MTL_BLAS_NAME(cgemm)(const char* transa, const char* transb, 
	    const int* m, const int* n, const int* k,
	    const std::complex<float>* alpha,  const std::complex<float> *da,  const int* lda,
	    const std::complex<float> *db, const int* ldb, const std::complex<float>* dbeta,
	    std::complex<float> *dc, const int* ldc);

void MTL_BLAS_NAME(zgemm)(const char* transa, const char* transb, 
	    const int* m, const int* n, const int* k,
	    const std::complex<double>* alpha,  const std::complex<double> *da,  const int* lda,
	    const std::complex<double> *db, const int* ldb, const std::complex<double>* dbeta,
	    std::complex<double> *dc, const int* ldc);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // MTL_HAS_BLAS

#endif // MTL_BLAS_INCLUDE
