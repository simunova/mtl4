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
#include <boost/timer.hpp>

#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/operation/cholesky.hpp>
#include <boost/numeric/mtl/operation/matrix_mult.hpp>
#include <boost/numeric/mtl/matrix/hessian_setup.hpp>
#include <boost/numeric/mtl/operation/assign_mode.hpp>
#include <boost/numeric/mtl/utility/papi.hpp>

using namespace mtl;
using namespace mtl::recursion; 
using namespace std;  



// Maximum time for a single measurement
// is 5 min
const double max_time= 300;

    // Bitmasks: 
    const unsigned long morton_mask= generate_mask<true, 0, row_major, 0>::value,
	morton_z_mask= generate_mask<false, 0, row_major, 0>::value,
	doppled_16_row_mask= generate_mask<true, 4, row_major, 0>::value,
	doppled_16_col_mask= generate_mask<true, 4, col_major, 0>::value,
	doppled_32_row_mask= generate_mask<true, 5, row_major, 0>::value,
	doppled_32_col_mask= generate_mask<true, 5, col_major, 0>::value,
	doppled_z_32_row_mask= generate_mask<false, 5, row_major, 0>::value,
	doppled_z_32_col_mask= generate_mask<false, 5, col_major, 0>::value,
	doppled_64_row_mask= generate_mask<true, 6, row_major, 0>::value,
	doppled_64_col_mask= generate_mask<true, 6, col_major, 0>::value,
	doppled_z_64_row_mask= generate_mask<false, 6, row_major, 0>::value,
	doppled_z_64_col_mask= generate_mask<false, 6, col_major, 0>::value,
	doppled_128_row_mask= generate_mask<true, 7, row_major, 0>::value,
	doppled_128_col_mask= generate_mask<true, 7, col_major, 0>::value,
	shark_32_row_mask= generate_mask<true, 5, row_major, 1>::value,
	shark_32_col_mask= generate_mask<true, 5, col_major, 1>::value,
	shark_z_32_row_mask= generate_mask<false, 5, row_major, 1>::value,
	shark_z_32_col_mask= generate_mask<false, 5, col_major, 1>::value,
	shark_64_row_mask= generate_mask<true, 6, row_major, 1>::value,
	shark_64_col_mask= generate_mask<true, 6, col_major, 1>::value,
	shark_z_64_row_mask= generate_mask<false, 6, row_major, 1>::value,
	shark_z_64_col_mask= generate_mask<false, 6, col_major, 1>::value;

typedef assign::plus_sum                            ama_t;

typedef recursion::bound_test_static<32>                    test32_t;
typedef recursion::bound_test_static<64>                    test64_t;

typedef gen_dense_mat_mat_mult_t<assign::plus_sum>  base_mult_t;
typedef gen_recursive_dense_mat_mat_mult_t<base_mult_t>     rec_mult_t;

typedef gen_tiling_22_dense_mat_mat_mult_t<assign::plus_sum>  tiling_22_base_mult_t;
typedef gen_tiling_44_dense_mat_mat_mult_t<assign::plus_sum>  tiling_44_base_mult_t;

utility::papi_t papi;
int l1i= papi.add_event("PAPI_L1_DCM");
int l2i= papi.add_event("PAPI_L2_DCM");
int tlbi= papi.add_event("PAPI_TLB_DM");


// ugly short cuts
typedef dense2D<double>                                       dr_t;
typedef dense2D<double, mat::parameters<col_major> >        dc_t;

#ifdef MTL_HAS_LAPACK
extern "C" {
  void dpotrf_(char* uplo, int*, double*, int*, int*);
}
#endif

struct dpotrf_t
{
    void operator()(dc_t& a)
    {
#ifdef MTL_HAS_LAPACK
	int size= a.num_rows(), info;
	dpotrf_("U", &size, &a[0][0], &size, &info);
#endif
    }
};



void print_time_and_mflops(double time, double size)
{
    // time and MFlops of single measure
    std::cout << time << ", " << 0.333333333333333333 * size * size * size / time / 1e6f << ", ";
}


// Matrices are only placeholder to provide the type
template <typename Visitor, typename Matrix>
void single_measure(Matrix&, Visitor visitor, unsigned size, std::vector<int>& enabled, int i)
{
    Matrix matrix(size, size);
    fill_matrix_for_cholesky(matrix);

    if (enabled[i]) {
	int reps= 0;
	boost::timer start;	
	papi.reset();
	for (; start.elapsed() < 5; reps++)
	    recursive_cholesky(matrix, visitor);
	papi.read();
	double time= start.elapsed() / double(reps);
	print_time_and_mflops(time, size);
	std::cout << papi[l1i]/reps << ", " << papi[l2i]/reps << ", " << papi[tlbi]/reps << ", ";
	if (time > max_time)
	    enabled[i]= 0;
    } else
	std::cout << ", , , , , ";
}


// Matrices are only placeholder to provide the type
template <typename Matrix, typename Functor>
void single_measure_dpotrf(Matrix&, Functor fcholesky, unsigned size, std::vector<int>& enabled, int i)
{
    Matrix matrix(size, size);
    fill_matrix_for_cholesky(matrix);

    if (enabled[i]) {
	int reps= 0;
	boost::timer start;	
	papi.reset();
	for (; start.elapsed() < 5; reps++)
	    fcholesky(matrix);
	papi.read();
	double time= start.elapsed() / double(reps);
	print_time_and_mflops(time, size);
	std::cout << papi[l1i]/reps << ", " << papi[l2i]/reps << ", " << papi[tlbi]/reps << ", ";
	if (time > max_time)
	    enabled[i]= 0;
    } else
	std::cout << ", , , , , ";
}


template <typename Visitor>
void measure(unsigned size, std::vector<int>& enabled)
{
    
    // The matrices in the following functions are only place holders, the real matrices are used in single_measure
    morton_dense<double,  morton_mask>             md(4, 4);
    morton_dense<double,  morton_z_mask>           mzd(4, 4);
    dense2D<double>                                dr(4, 4);
    dense2D<double, mat::parameters<col_major> > dc(4, 4);
    morton_dense<double,  doppled_32_row_mask>     d32r(4, 4);
    morton_dense<double,  doppled_64_row_mask>     d64r(4, 4);

    std::cout << size << ", ";

    Visitor visitor;
    single_measure(md, visitor, size, enabled, 0);
    single_measure(mzd, visitor, size, enabled, 1);
    single_measure(dr, visitor, size, enabled, 2);
    single_measure(dc, visitor, size, enabled, 3);
    single_measure(d32r, visitor, size, enabled, 4);
    single_measure(d64r, visitor, size, enabled, 5);

    single_measure_dpotrf(dc, dpotrf_t(), size, enabled, 6);
    
    std::cout << "0\n"; // to not finish with comma
    std::cout.flush();
}


template <typename Measure>
void series(unsigned steps, unsigned max_size, Measure measure, const string& comment)
{
    std::cout << "# " << comment << '\n';
    std::cout << "# Gnu-Format size, time, MFlops, ...\n"; 
    std::cout << "# N-order, Z-order, row, col, 32 row, 64 row\n"; std::cout.flush();

    std::vector<int> enabled(16, 1);
    for (unsigned i= steps; i <= max_size; i+= steps)
	measure(i, enabled);
}

int main(int argc, char* argv[])
{
    papi.start();

    std::vector<std::string> scenarii;
    scenarii.push_back(string("Comparing canonical implementation with different matrix types"));
    scenarii.push_back(string("Comparing iterator implementation with different matrix types"));
    scenarii.push_back(string("Comparing prev. using fast Schur update (2x2 tiling) with different matrix types"));
    scenarii.push_back(string("Comparing prev. using fast Schur update (4x4 tiling) with different matrix types"));

    using std::cout;
    if (argc < 4) {
	cerr << "usage: recursive_mult_timing <scenario> <steps> <max_size>\nScenarii:\n"; 
	for (unsigned i= 0; i < scenarii.size(); i++)
	    cout << i << ": " << scenarii[i] << "\n";
	exit(1);
    }
    unsigned int scenario= atoi(argv[1]), steps= atoi(argv[2]), max_size= atoi(argv[3]), size= 32; 

    typedef with_bracket::recursive_cholesky_base_visitor_t     bracket_t;
    typedef with_iterator::recursive_cholesky_base_visitor_t    iterator_t;

    typedef detail::mult_schur_update_t<gen_tiling_22_dense_mat_mat_mult_t<assign::minus_sum> > schur_update_22_t;
    typedef recursive_cholesky_visitor_t<recursion::bound_test_static<64>, with_iterator::cholesky_base_t, with_iterator::tri_solve_base_t, 
	                                 with_iterator::tri_schur_base_t, schur_update_22_t>   
	tiling_22_t;

    typedef detail::mult_schur_update_t<gen_tiling_44_dense_mat_mat_mult_t<assign::minus_sum> > schur_update_44_t;
    typedef recursive_cholesky_visitor_t<recursion::bound_test_static<64>, with_iterator::cholesky_base_t, with_iterator::tri_solve_base_t, 
	                                 with_iterator::tri_schur_base_t, schur_update_44_t>   
	tiling_44_t;


    switch (scenario) {
      case 0: 	series(steps, max_size, measure<bracket_t>, scenarii[0]); break;
      case 1: 	series(steps, max_size, measure<iterator_t>, scenarii[1]); break;
      case 2: 	series(steps, max_size, measure<tiling_22_t>, scenarii[2]); break;
      case 3: 	series(steps, max_size, measure<tiling_44_t>, scenarii[3]); break;
    }

    return 0; 

}
 


// Compile command:
// g++4 cholesky_timing.cpp -o cholesky_timing -O3 -DNDEBUG -ffast-math  -mcpu=opteron -mtune=opteron -msse2 -mfpmath=sse  -I${MTL_BOOST_ROOT} -I${BOOST_ROOT} -I/usr/local/include -L/usr/local/lib -lpapi -DMTL_HAS_PAPI -DMTL_HAS_BLAS -DMTL_HAS_LAPACK -L/u/htor/projekte/mathlibs/acml-2-6-0-gnu-64bit/gnu64/lib -lacml -L/usr/lib/gcc/x86_64-redhat-linux/3.4.3 -lg2c
