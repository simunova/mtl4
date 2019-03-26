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

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include <boost/timer.hpp>

#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/operation/matrix_mult.hpp>
#include <boost/numeric/mtl/matrix/hessian_setup.hpp>
#include <boost/numeric/mtl/operation/assign_mode.hpp>
#include <boost/numeric/mtl/recursion/predefined_masks.hpp>
#include <boost/numeric/mtl/utility/papi.hpp>

using namespace mtl;
using namespace mtl::recursion; 
using namespace std;  

// Maximum time for a single measurement
// is 20 min
const double max_time= 900;

typedef assign::plus_sum                            ama_t;

typedef recursion::bound_test_static<32>                    test32_t;
typedef recursion::bound_test_static<64>                    test64_t;

typedef gen_dense_mat_mat_mult_t<assign::plus_sum>  base_mult_t;
typedef gen_recursive_dense_mat_mat_mult_t<base_mult_t>     rec_mult_t;

typedef gen_tiling_22_dense_mat_mat_mult_t<assign::plus_sum>  tiling_22_base_mult_t;
typedef gen_tiling_44_dense_mat_mat_mult_t<assign::plus_sum>  tiling_44_base_mult_t;

// ugly short cuts
typedef dense2D<double>                                       dr_t;
typedef dense2D<double, mat::parameters<col_major> >        dc_t;

utility::papi_t papi;
int l1i= papi.add_event("PAPI_L1_DCM");
int l2i= papi.add_event("PAPI_L2_DCM");
int tlbi= papi.add_event("PAPI_TLB_DM");

#ifdef MTL_USE_BLAS
extern "C" {
void dgemm_(const char* transa, const char* transb, 
	    const int* m, const int* n, const int* k,
	    const double* alpha,  const double *da,  const int* lda,
	    const double *db, const int* ldb, const double* dbeta,
	    double *dc, const int* ldc);
}
#endif

struct dgemm_t
{
    void operator()(const dc_t& a, const dc_t& b, dc_t& c)
    {
#ifdef MTL_USE_BLAS
	int size= a.num_rows();
	double alpha= 1.0, beta= 0.0;
	dgemm_("N", "N", &size, &size, &size, &alpha, 
	       const_cast<double*>(&a[0][0]), &size, const_cast<double*>(&b[0][0]), 
	       &size, &beta, &c[0][0], &size);
#endif
    }
};

struct dgemm_add_t
{
    void operator()(const dc_t& a, const dc_t& b, dc_t& c)
    {
	int size= a.num_rows();
	double alpha= 1.0, beta= 1.0;
	dgemm_("N", "N", &size, &size, &size, &alpha, 
	       const_cast<double*>(&a[0][0]), &size, const_cast<double*>(&b[0][0]), 
	       &size, &beta, &c[0][0], &size);

    }
};



void print_time_and_mflops(double time, double size)
{
    // time and MFlops of single measure
    std::cout << time << ", " << 2.0 * size * size * size / time / 1e6f << ", ";
}


// Matrices are only placeholder to provide the type
template <typename MatrixA, typename MatrixB, typename MatrixC, typename Mult>
void single_measure(MatrixA&, MatrixB&, MatrixC&, Mult mult, unsigned size, std::vector<int>& enabled, int i)
{
    MatrixA a(size, size);
    MatrixB b(size, size);
    MatrixC c(size, size);
    hessian_setup(a, 1.0);
    hessian_setup(b, 1.0);

    if (enabled[i]) {
	int reps= 0;
	boost::timer start;	
	papi.reset();
	for (; start.elapsed() < 5; reps++)
	    mult(a, b, c);
	papi.read();
	double time= start.elapsed() / double(reps);

	print_time_and_mflops(time, a.num_rows());
	std::cout << papi[l1i]/reps << ", " << papi[l2i]/reps << ", " << papi[tlbi]/reps << ", ";
	if (time > max_time)
	    enabled[i]= 0;
    } else
	std::cout << ", , , , , ";
}

// The matrices in the following functions are only place holders, the real matrices are used in single_measure
void measure_morton_order(unsigned size, std::vector<int>& enabled)
{
    morton_dense<double,  morton_mask>             mda(4, 4), mdb(4, 4), mdc(4, 4);
    morton_dense<double,  morton_z_mask>           mzda(4, 4), mzdb(4, 4), mzdc(4, 4);
    dc_t                                           dc(4, 4);
    
    std::cout << size << ", ";
    rec_mult_t  mult;
    single_measure(mda, mdb, mdc, mult, size, enabled, 0);
    single_measure(mzda, mzdb, mzdc, mult, size, enabled, 1);
    single_measure(mda, mzdb, mdc, mult, size, enabled, 2);
    single_measure(dc, dc, dc, dgemm_t(), size, enabled, 3);
    
    std::cout << "0\n"; // to not finish with comma
    std::cout.flush();
}


void measure_cast(unsigned size, std::vector<int>& enabled)
{
    morton_dense<double,  morton_mask>             mda(4, 4), mdb(4, 4), mdc(4, 4);
    morton_dense<double,  doppled_32_row_mask>     d32ra(4, 4), d32rb(4, 4), d32rc(4, 4);
    morton_dense<double,  doppled_64_row_mask>     d64ra(4, 4), d64rb(4, 4), d64rc(4, 4);
    morton_dense<double,  doppled_64_col_mask>     d64ca(4, 4), d64cb(4, 4), d64cc(4, 4); 
    
    rec_mult_t  mult;
    std::cout << size << ", ";
    single_measure(mda, mdb, mdc, mult, size, enabled, 0);
    single_measure(d32ra, d32rb, d32rc, mult, size, enabled, 1);
    single_measure(d64ra, d64rb, d64rc, mult, size, enabled, 2);
    single_measure(d64ca, d64cb, d64cc, mult, size, enabled, 3);

    dc_t                                           dc(4, 4);
    single_measure(dc, dc, dc, dgemm_t(), size, enabled, 4);
 
    std::cout << "0\n";  std::cout.flush();
}


void measure_with_unroll(unsigned size, std::vector<int>& enabled)
{
    morton_dense<double,  doppled_32_row_mask>     d32r(4, 4);
    morton_dense<double,  doppled_32_col_mask>     d32c(4, 4);

    std::cout << size << ", ";

    gen_recursive_dense_mat_mat_mult_t<base_mult_t, test32_t>           mult;
    gen_recursive_dense_mat_mat_mult_t<tiling_22_base_mult_t, test32_t> mult_22;
    gen_recursive_dense_mat_mat_mult_t<tiling_44_base_mult_t, test32_t> mult_44;

    single_measure(d32r, d32r, d32r, mult, size, enabled, 0);
    single_measure(d32r, d32r, d32r, mult_22, size, enabled, 1);
    single_measure(d32r, d32r, d32r, mult_44, size, enabled, 2);

    single_measure(d32c, d32c, d32c, mult, size, enabled, 3);
    single_measure(d32c, d32c, d32c, mult_22, size, enabled, 4);
    single_measure(d32c, d32c, d32c, mult_44, size, enabled, 5);
    dc_t                                           dc(4, 4);
    single_measure(dc, dc, dc, dgemm_t(), size, enabled, 6);
 
    std::cout << "0\n";  std::cout.flush();
}

void measure_base_size(unsigned size, std::vector<int>& enabled)
{
    morton_dense<double,  doppled_16_row_mask>     d16r(4, 4);
    morton_dense<double,  doppled_32_row_mask>     d32r(4, 4);
    morton_dense<double,  doppled_64_row_mask>     d64r(4, 4);
    morton_dense<double,  doppled_128_col_mask>    d128r(4, 4);
    
    std::cout << size << ", ";

    gen_recursive_dense_mat_mat_mult_t<tiling_22_base_mult_t, recursion::bound_test_static<16> > mult16;
    single_measure(d16r, d16r, d16r, mult16, size, enabled, 0);

    gen_recursive_dense_mat_mat_mult_t<tiling_22_base_mult_t, recursion::bound_test_static<32> > mult32;
    single_measure(d32r, d32r, d32r, mult32, size, enabled, 1);

    gen_recursive_dense_mat_mat_mult_t<tiling_22_base_mult_t, recursion::bound_test_static<64> > mult64;
    single_measure(d64r, d64r, d64r, mult64, size, enabled, 2);

    gen_recursive_dense_mat_mat_mult_t<tiling_22_base_mult_t, recursion::bound_test_static<128> > mult128;
    single_measure(d128r, d128r, d128r, mult128, size, enabled, 3);

    dc_t                                           dc(4, 4);
    single_measure(dc, dc, dc, dgemm_t(), size, enabled, 4);
 
    std::cout << "0\n";  std::cout.flush();
}

template <typename Matrix, typename MatrixB> 
void measure_unrolling(unsigned size, std::vector<int>& enabled, Matrix& matrix, MatrixB& matrixb)
{
    std::cout << size << ", ";
 
    gen_recursive_dense_mat_mat_mult_t<base_mult_t>           mult;
    gen_recursive_dense_mat_mat_mult_t<tiling_22_base_mult_t> mult_22;
    gen_recursive_dense_mat_mat_mult_t<tiling_44_base_mult_t> mult_44;

    typedef gen_tiling_dense_mat_mat_mult_t<2, 2, ama_t>  tiling_m22_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m22_base_mult_t> mult_m22;
    
    typedef gen_tiling_dense_mat_mat_mult_t<2, 4, ama_t>  tiling_m24_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m24_base_mult_t> mult_m24;

    typedef gen_tiling_dense_mat_mat_mult_t<4, 2, ama_t>  tiling_m42_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m42_base_mult_t> mult_m42;

    typedef gen_tiling_dense_mat_mat_mult_t<3, 5, ama_t>  tiling_m35_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m35_base_mult_t> mult_m35;

    typedef gen_tiling_dense_mat_mat_mult_t<4, 4, ama_t>  tiling_m44_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m44_base_mult_t> mult_m44;


    single_measure(matrix, matrixb, matrix, mult, size, enabled, 0);
    single_measure(matrix, matrixb, matrix, mult_22, size, enabled, 1);
    single_measure(matrix, matrixb, matrix, mult_44, size, enabled, 2);

    single_measure(matrix, matrixb, matrix, mult_m22, size, enabled, 3);
    single_measure(matrix, matrixb, matrix, mult_m24, size, enabled, 4);
    single_measure(matrix, matrixb, matrix, mult_m42, size, enabled, 5);
    single_measure(matrix, matrixb, matrix, mult_m35, size, enabled, 6);
    single_measure(matrix, matrixb, matrix, mult_m44, size, enabled, 7);
 
    dc_t                                           dc(4, 4);
    single_measure(dc, dc, dc, dgemm_t(), size, enabled, 8);

    gen_recursive_dense_mat_mat_mult_t<dgemm_add_t> mult_blas;
    single_measure(dc, dc, dc, mult_blas, size, enabled, 9);


    std::cout << "0\n";  std::cout.flush();
}

void measure_unrolling_hybrid(unsigned size, std::vector<int>& enabled)
{
    morton_dense<double,  doppled_64_row_mask>     d64r(4, 4);
    morton_dense<double,  doppled_64_col_mask>     d64c(4, 4);
    measure_unrolling(size, enabled, d64r, d64c);
}

void measure_unrolling_dense(unsigned size, std::vector<int>& enabled)
{
    dense2D<double> dense(4, 4);
    dense2D<double, mat::parameters<col_major> >    b(4, 4);
    measure_unrolling(size, enabled, dense, b);
}


void measure_orientation(unsigned size, std::vector<int>& enabled)
{
    std::cout << size << ", ";
 
    morton_dense<double,  doppled_64_row_mask>     d64r(4, 4);
    morton_dense<double,  doppled_64_col_mask>     d64c(4, 4);

    typedef gen_tiling_dense_mat_mat_mult_t<4, 2, ama_t>  tiling_m42_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m42_base_mult_t> mult_m42;

    typedef gen_tiling_dense_mat_mat_mult_t<4, 4, ama_t>  tiling_m44_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m44_base_mult_t> mult_m44;
    
    single_measure(d64r, d64r, d64r, mult_m42, size, enabled, 0);
    single_measure(d64c, d64c, d64c, mult_m42, size, enabled, 1);
    single_measure(d64r, d64c, d64r, mult_m42, size, enabled, 2);
    single_measure(d64c, d64r, d64r, mult_m42, size, enabled, 3);
    
    single_measure(d64r, d64r, d64r, mult_m44, size, enabled, 4);
    single_measure(d64c, d64c, d64c, mult_m44, size, enabled, 5);
    single_measure(d64r, d64c, d64r, mult_m44, size, enabled, 6);
    single_measure(d64c, d64r, d64r, mult_m44, size, enabled, 7);
 
    dc_t                                           dc(4, 4);
    single_measure(dc, dc, dc, dgemm_t(), size, enabled, 8);

    std::cout << "0\n";  std::cout.flush();
}


void measure_unrolling_32(unsigned size, std::vector<int>& enabled)
{
    std::cout << size << ", ";
 
    gen_recursive_dense_mat_mat_mult_t<base_mult_t, test32_t>           mult;
    gen_recursive_dense_mat_mat_mult_t<tiling_22_base_mult_t, test32_t> mult_22;
    gen_recursive_dense_mat_mat_mult_t<tiling_44_base_mult_t, test32_t> mult_44;

    typedef gen_tiling_dense_mat_mat_mult_t<2, 2, ama_t>  tiling_m22_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m22_base_mult_t, test32_t> mult_m22;
    
    typedef gen_tiling_dense_mat_mat_mult_t<2, 4, ama_t>  tiling_m24_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m24_base_mult_t, test32_t> mult_m24;

    typedef gen_tiling_dense_mat_mat_mult_t<4, 2, ama_t>  tiling_m42_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m42_base_mult_t, test32_t> mult_m42;

    typedef gen_tiling_dense_mat_mat_mult_t<3, 5, ama_t>  tiling_m35_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m35_base_mult_t, test32_t> mult_m35;

    typedef gen_tiling_dense_mat_mat_mult_t<4, 4, ama_t>  tiling_m44_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m44_base_mult_t, test32_t> mult_m44;

    
    morton_dense<double,  doppled_32_row_mask>     d32r(4, 4);
    morton_dense<double,  doppled_32_col_mask>     d32c(4, 4);


    single_measure(d32r, d32c, d32r, mult, size, enabled, 0);
    single_measure(d32r, d32c, d32r, mult_22, size, enabled, 1);
    single_measure(d32r, d32c, d32r, mult_44, size, enabled, 2);

    single_measure(d32r, d32c, d32r, mult_m22, size, enabled, 3);
    single_measure(d32r, d32c, d32r, mult_m24, size, enabled, 4);
    single_measure(d32r, d32c, d32r, mult_m42, size, enabled, 5);
    single_measure(d32r, d32c, d32r, mult_m35, size, enabled, 6);
    single_measure(d32r, d32c, d32r, mult_m44, size, enabled, 7);
 
    dc_t                                           dc(4, 4);
    single_measure(dc, dc, dc, dgemm_t(), size, enabled, 8);

    std::cout << "0\n";  std::cout.flush();
}



void measure_hetero_value(unsigned size, std::vector<int>& enabled)
{
    using std::complex;
    typedef gen_tiling_dense_mat_mat_mult_t<4, 4, ama_t>  tiling_m44_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m44_base_mult_t> mult;

    dc_t                                           dc(4, 4);
    dr_t                                           dr(4, 4);
    dense2D<float, mat::parameters<row_major> >  fr(4, 4);
    dense2D<float, mat::parameters<col_major> >  fc(4, 4);
    dense2D<complex<float>, mat::parameters<col_major> >   cc(4, 4);
    dense2D<complex<double>, mat::parameters<col_major> >  zc(4, 4);
    dense2D<complex<double>, mat::parameters<row_major> >  zr(4, 4);

    morton_dense<double,  doppled_64_row_mask>     d64r(4, 4);
    morton_dense<double,  doppled_64_col_mask>     d64c(4, 4);
    
    morton_dense<float,  doppled_64_row_mask>     f64r(4, 4);
    morton_dense<float,  doppled_64_col_mask>     f64c(4, 4);

    morton_dense<complex<float>,  doppled_64_row_mask>     c64r(4, 4);
    morton_dense<complex<float>,  doppled_64_col_mask>     c64c(4, 4);

    morton_dense<complex<double>,  doppled_64_row_mask>     z64r(4, 4);
    morton_dense<complex<double>,  doppled_64_col_mask>     z64c(4, 4);

    std::cout << size << ", ";

    single_measure(dr, fc, dr, mult, size, enabled, 0);
    single_measure(fr, fc, dr, mult, size, enabled, 1);
   
    //single_measure(fr, zc, zr, mult, size, enabled, 2);
    //single_measure(dr, cc, zr, mult, size, enabled, 3);
    
    single_measure(d64r, f64c, d64r, mult, size, enabled, 4);
    single_measure(d64r, z64c, z64r, mult, size, enabled, 5);
    std::cout << "0\n";  std::cout.flush();
}


void measure_hetero_layout(unsigned size, std::vector<int>& enabled)
{
    using std::complex;
    typedef gen_tiling_dense_mat_mat_mult_t<4, 4, ama_t>  tiling_m44_base_mult_t;
    gen_recursive_dense_mat_mat_mult_t<tiling_m44_base_mult_t> mult;

    dc_t                                           dc(4, 4);
    dr_t                                           dr(4, 4);
    dense2D<float, mat::parameters<row_major> >  fr(4, 4);
    dense2D<float, mat::parameters<col_major> >  fc(4, 4);
    dense2D<complex<float>, mat::parameters<col_major> >   cc(4, 4);
    dense2D<complex<double>, mat::parameters<col_major> >  zc(4, 4);
    dense2D<complex<double>, mat::parameters<row_major> >  zr(4, 4);

    morton_dense<double,  doppled_64_row_mask>     d64r(4, 4);
    morton_dense<double,  doppled_64_col_mask>     d64c(4, 4);
    
    morton_dense<float,  doppled_64_row_mask>     f64r(4, 4);
    morton_dense<float,  doppled_64_col_mask>     f64c(4, 4);

    morton_dense<complex<float>,  doppled_64_row_mask>     c64r(4, 4);
    morton_dense<complex<float>,  doppled_64_col_mask>     c64c(4, 4);

    morton_dense<complex<double>,  doppled_64_row_mask>     z64r(4, 4);
    morton_dense<complex<double>,  doppled_64_col_mask>     z64c(4, 4);

    std::cout << size << ", ";

    single_measure(dr, f64c, dr, mult, size, enabled, 0);
    single_measure(fr, f64c, dr, mult, size, enabled, 1);
   
    //single_measure(f64r, zc, zr, mult, size, enabled, 2);
    single_measure(d64r, zc, zr, mult, size, enabled, 3);
    
    single_measure(d64r, f64c, dr, mult, size, enabled, 4);
    //single_measure(f64r, z64c, zr, mult, size, enabled, 5);
    std::cout << "0\n";  std::cout.flush();
}




template <typename Measure>
void series(unsigned steps, unsigned max_size, Measure measure, const string& comment)
{
    std::cout << "# " << comment << '\n';
    std::cout << "# Gnu-Format size, (time, MFlops, L1, L2, TLB, time, MFlops, )*\n"; std::cout.flush();

    std::vector<int> enabled(16, 1);
    for (unsigned i= steps; i <= max_size; i+= steps)
	measure(i, enabled);
}



int main(int argc, char* argv[])
{
    papi.start();

    std::vector<std::string> scenarii;
    scenarii.push_back(string("Comparing Z-, N-order and mixed with recursive multiplication"));
    scenarii.push_back(string("Comparing base case cast (64) for Z-order, hybrid 32, hybrid 64 row and col-major"));
    scenarii.push_back(string("Using unrolled mult on hybrid row- and column-major matrices"));
    scenarii.push_back(string("Comparing base case sizes for corresponding hybrid row-major matrices"));
    scenarii.push_back(string("Comparing different unrolling for hybrid row-major matrices"));
    scenarii.push_back(string("Comparing different unrolling for row-major dense matrices"));
    scenarii.push_back(string("Comparing different orientations for hybrid row-major matrices"));
    scenarii.push_back(string("Comparing different unrolling for hybrid 32 row-major times col-major matrices"));
    scenarii.push_back(string("Multiplying matrices with different value types"));
    scenarii.push_back(string("Multiplying matrices with different matrix layouts"));

    using std::cout;
    if (argc < 4) {
	cerr << "usage: recursive_mult_timing <scenario> <steps> <max_size>\nScenarii:\n"; 
	for (unsigned i= 0; i < scenarii.size(); i++)
	    cout << i << ": " << scenarii[i] << "\n";
	exit(1);
    }
    unsigned int scenario= atoi(argv[1]), steps= atoi(argv[2]), max_size= atoi(argv[3]), size= 32; 

    switch (scenario) {
      case 0: 	series(steps, max_size, measure_morton_order, scenarii[0]); break;
      case 1: 	series(steps, max_size, measure_cast, scenarii[1]); break;
      case 2: 	series(steps, max_size, measure_with_unroll, scenarii[2]); break;
      case 3: 	series(steps, max_size, measure_base_size, scenarii[3]); break;
      case 4: 	series(steps, max_size, measure_unrolling_hybrid, scenarii[4]); break;
      case 5: 	series(steps, max_size, measure_unrolling_dense, scenarii[5]); break;
      case 6: 	series(steps, max_size, measure_orientation, scenarii[6]); break;
      case 7: 	series(steps, max_size, measure_unrolling_32, scenarii[7]); break;
      case 8: 	series(steps, max_size, measure_hetero_value, scenarii[8]); break;
      case 9: 	series(steps, max_size, measure_hetero_layout, scenarii[9]); break;
    }

    return 0; 

}
 





#if 0

// scheiss Kommandos

// g++4 matrix_product_timing.cpp  -o matrix_product_timing -O3 -DNDEBUG -ffast-math  -mcpu=opteron -mtune=opteron -msse2 -mfpmath=sse -I${MTL_BOOST_ROOT} -I${BOOST_ROOT} -I/usr/local/include -L/usr/local/lib -lpapi -DMTL_HAS_PAPI -DMTL_HAS_BLAS  -L/u/htor/projekte/mathlibs/goto-blas -lgoto_opteron-64 -lpthread xerbla.o


// g++4 matrix_product_timing.cpp  -o matrix_product_timing -O3 -DNDEBUG -ffast-math  -mcpu=opteron -mtune=opteron -msse2 -mfpmath=sse -I${MTL_BOOST_ROOT} -I${BOOST_ROOT} -I/usr/local/include -L/usr/local/lib -lpapi -DMTL_HAS_PAPI -DMTL_HAS_BLAS  -L/u/htor/projekte/mathlibs/acml-2-6-0-gnu-64bit/gnu64/lib -lacml -L/usr/lib/gcc/x86_64-redhat-linux/3.4.3 -lg2c

// g++4 matrix_product_timing.cpp -o matrix_product_timing -O3 -DNDEBUG -ffast-math -mcpu=opteron -mtune=opteron -msse2 -mfpmath=sse -I${MTL_BOOST_ROOT} -I${BOOST_ROOT} -I/usr/local/include -L/usr/local/lib -lpapi -DMTL_HAS_PAPI -DMTL_HAS_BLAS -L/san/atlas/lib -lf77blas -latlas -L/usr/lib/gcc/x86_64-redhat-linux/3.4.3 -lg2c

#endif
