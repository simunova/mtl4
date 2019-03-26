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
#include <vector>
#include <boost/test/minimal.hpp>
#include <boost/timer.hpp>
#include <boost/static_assert.hpp>

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/hessian_setup.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/operation/assign_mode.hpp>

using namespace std;
using namespace mtl;

template <unsigned Outer, unsigned OuterMax, unsigned Inner, unsigned InnerMax>
struct double_unroll
{
    static unsigned const outer= Outer - 1, next_outer= Outer,
            	          inner= Inner - 1, next_inner= Inner + 1;
    void operator() ()
    {
	cout << outer << " : " << inner << "\n";
	double_unroll<next_outer, OuterMax, next_inner, InnerMax>() ();
    }
};


template <unsigned Outer, unsigned OuterMax, unsigned InnerMax>
struct double_unroll<Outer, OuterMax, InnerMax, InnerMax>
{
    static unsigned const outer= Outer - 1, next_outer= Outer + 1,
            	          inner= InnerMax - 1, next_inner= 1;
    void operator() ()
    {
	cout << outer << " : " << inner << "\n";
	double_unroll<next_outer, OuterMax, next_inner, InnerMax>() ();
    }
};


template <unsigned OuterMax, unsigned InnerMax>
struct double_unroll<OuterMax, OuterMax, InnerMax, InnerMax>
{
    static unsigned const outer= OuterMax - 1,
            	          inner= InnerMax - 1;
    void operator() ()
    {
	cout << outer << " : " << inner << "\n";
    }
};

template <unsigned Outer, unsigned OuterMax, unsigned Middle, unsigned MiddleMax,
	  unsigned Inner, unsigned InnerMax>
struct triple_unroll
{
    static unsigned const outer= Outer - 1, next_outer= Outer,
	                  middle= Middle - 1, next_middle= Middle,
            	          inner= Inner - 1, next_inner= Inner + 1;
    void operator() ()
    {
	cout << outer << " : " << middle << " : " << inner << "\n";
	triple_unroll<next_outer, OuterMax, next_middle, MiddleMax, next_inner, InnerMax>() ();
    }
};


template <unsigned Outer, unsigned OuterMax, unsigned Middle, unsigned MiddleMax, unsigned InnerMax>
struct triple_unroll<Outer, OuterMax, Middle, MiddleMax, InnerMax, InnerMax>
{
    static unsigned const outer= Outer - 1, next_outer= Outer,
	                  middle= Middle - 1, next_middle= Middle + 1,
            	          inner= InnerMax - 1, next_inner= 1;
    void operator() ()
    {
	cout << outer << " : " << middle << " : " << inner << "\n";
	triple_unroll<next_outer, OuterMax, next_middle, MiddleMax, next_inner, InnerMax>() ();
    }
};


template <unsigned Outer, unsigned OuterMax, unsigned MiddleMax, unsigned InnerMax>
struct triple_unroll<Outer, OuterMax, MiddleMax, MiddleMax, InnerMax, InnerMax>
{
    static unsigned const outer= Outer - 1, next_outer= Outer + 1,
	                  middle= MiddleMax - 1, next_middle= 1,
            	          inner= InnerMax - 1, next_inner= 1;
    void operator() ()
    {
	cout << outer << " : " << middle << " : " << inner << "\n";
	triple_unroll<next_outer, OuterMax, next_middle, MiddleMax, next_inner, InnerMax>() ();
    }
};


template <unsigned OuterMax, unsigned MiddleMax, unsigned InnerMax>
struct triple_unroll<OuterMax, OuterMax, MiddleMax, MiddleMax, InnerMax, InnerMax>
{
    static unsigned const outer= OuterMax - 1,
	                  middle= MiddleMax - 1,
            	          inner= InnerMax - 1;
    void operator() ()
    {
	cout << outer << " : " << middle << " : " << inner << "\n";
    }
};


template <unsigned Outer, unsigned OuterMax, unsigned Middle, unsigned MiddleMax,
	  unsigned Inner, unsigned InnerMax>
struct triple_wrapper : public triple_unroll<Outer, OuterMax, Middle, MiddleMax, Inner, InnerMax>
{
    typedef triple_unroll<Outer, OuterMax, Middle, MiddleMax, Inner, InnerMax> base;
    typedef triple_wrapper<base::next_outer, OuterMax, base::next_middle, MiddleMax, base::next_inner, InnerMax> next_t;

    void operator() ()
    {
	cout << this->outer << " : " << this->middle << " : " << this->inner << "\n";
	next_t() ();
    }  
};


template <unsigned OuterMax, unsigned MiddleMax, unsigned InnerMax>
struct triple_wrapper<OuterMax, OuterMax, MiddleMax, MiddleMax, InnerMax, InnerMax>
    : public triple_unroll<OuterMax, OuterMax, MiddleMax, MiddleMax, InnerMax, InnerMax>
{
    void operator() ()
    {
	cout << this->outer << " : " << this->middle << " : " << this->inner << "\n";
    }  
};



void print_time_and_mflops(double time, double size)
{ 
    // std::cout << "    takes " << time << "s = " << 2.0 * size * size * size / time / 1e6f << "MFlops\n";
    std::cout << size << ", " << time << ", " << 2.0 * size * size * size / time / 1e6f << "\n";
    std::cout.flush();
}


// Matrices are only placeholder to provide the type
template <typename MatrixA, typename MatrixB, typename MatrixC, typename Mult>
double time_measure(MatrixA&, MatrixB&, MatrixC&, Mult mult, unsigned size)
{
    MatrixA a(size, size);
    MatrixB b(size, size);
    MatrixC c(size, size);

    hessian_setup(a, 1.0);
    hessian_setup(b, 2.0); 

    // repeat multiplication if it is less than a second (until it is a second)
    int i; boost::timer start1;
    for (i= 0; start1.elapsed() < 1.0; i++)
	mult(a, b, c);
    double elapsed= start1.elapsed() / double(i);
    print_time_and_mflops(elapsed, size);
    check_hessian_matrix_product(c, size);
    return elapsed;
}
 
template <typename MatrixA, typename MatrixB, typename MatrixC, typename Mult>
void time_series(MatrixA& a, MatrixB& b, MatrixC& c, Mult mult, const char* name, unsigned steps, unsigned max_size)
{
    // Maximal time per measurement 20 min
    double max_time= 1200.0;

    std::cout << "\n# " << name << ":\n";
    std::cout << "# Gnu-Format size, time, MFlops\n";
    std::cout.flush();

    for (unsigned i= steps; i <= max_size; i+= steps) {
	double elapsed= time_measure(a, b, c, mult, i);
	if (elapsed > max_time) break;
    }
}


void mult_simple(const dense2D<double>& a, const dense2D<double>& b, dense2D<double>& c)
{
    for (unsigned i= 0; i < c.num_rows(); i++)
	for (unsigned k= 0; k < c.num_cols(); k++) {
	    double tmp= 0.0;
	    for (unsigned j= 0; j < b.num_cols(); j++)
		tmp+= a[i][j] * b[j][k];
	    c[i][k]= tmp;
	}
}

void mult_simple_p(dense2D<double>& a, dense2D<double>& b, dense2D<double>& c)
{
    for (unsigned i= 0; i < c.num_rows(); i++)
	for (unsigned k= 0; k < c.num_cols(); k++) {
	    double tmp= 0.0;
	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()], *begin_b= &b[0][k];
	    int ld= b.num_rows();
	    for (; begin_a != end_a; ++begin_a, begin_b+= ld)
		tmp+= *begin_a * *begin_b;
	    c[i][k]= tmp;
	}
}

typedef dense2D<double, mat::parameters<col_major> >   cm_type;
typedef dense2D<double, mat::parameters<row_major> >   rm_type;

void mult_simple_pt(dense2D<double>& a, cm_type& b, dense2D<double>& c)
{
    for (unsigned i= 0; i < c.num_rows(); i++)
	for (unsigned k= 0; k < c.num_cols(); k++) {
	    double tmp= 0.0;
	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()], *begin_b= &b[0][k];
	    for (; begin_a != end_a; ++begin_a, ++begin_b)
		tmp+= *begin_a * *begin_b;
	    c[i][k]= tmp;
	}
}


void mult_simple_ptu(dense2D<double>& a, cm_type& b, dense2D<double>& c)
{
    for (unsigned i= 0; i < c.num_rows(); i++)
	for (unsigned k= 0; k < c.num_cols(); k+=2) {
	    int ld= b.num_rows();
	    double tmp0= 0.0, tmp1= 0.0;

	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()];
	    double *begin_b= &b[0][k];
	    for (; begin_a != end_a; ++begin_a, ++begin_b) {
		tmp0+= *begin_a * *begin_b;
		tmp1+= *begin_a * *(begin_b+ld);
	    }
	    c[i][k]= tmp0; c[i][k+1]= tmp1;
	}
}

void mult_simple_ptu4(dense2D<double>& a, cm_type& b, dense2D<double>& c)
{
    for (unsigned i= 0; i < c.num_rows(); i++)
	for (unsigned k= 0; k < c.num_cols(); k+=4) {
	    int ld1= b.num_rows(), ld2= 2*ld1, ld3=3*ld1;
	    double tmp0= 0.0, tmp1= 0.0, tmp2= 0.0, tmp3= 0.0;

	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()];
	    double *begin_b= &b[0][k];
	    for (; begin_a != end_a; ++begin_a, ++begin_b) {
		tmp0+= *begin_a * *begin_b;
		tmp1+= *begin_a * *(begin_b+ld1);
		tmp2+= *begin_a * *(begin_b+ld2);
		tmp3+= *begin_a * *(begin_b+ld3);
	    }
	    c[i][k]= tmp0; c[i][k+1]= tmp1;
	    c[i][k+2]= tmp2; c[i][k+3]= tmp3;
	}
}

// C must be square!
void mult_simple_ptu22(dense2D<double>& a, cm_type& b, dense2D<double>& c)
{
    for (unsigned i= 0; i < c.num_rows(); i+=2)
	for (unsigned k= 0; k < c.num_cols(); k+=2) {
	    int ld= b.num_rows();
	    double tmp00= 0.0, tmp01= 0.0, tmp10= 0.0, tmp11= 0.0;

	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()];
	    double *begin_b= &b[0][k];
	    for (; begin_a != end_a; ++begin_a, ++begin_b) {
		tmp00+= *begin_a * *begin_b;
		tmp01+= *begin_a * *(begin_b+ld);
		tmp10+= *(begin_a+ld) * *begin_b;
		tmp11+= *(begin_a+ld) * *(begin_b+ld);
	    }
	    c[i][k]= tmp00; c[i][k+1]= tmp01;
	    c[i+1][k]= tmp10; c[i+1][k+1]= tmp11;
	}
}

// C must have even dimensions
template <typename MatrixA, typename MatrixB, typename MatrixC>
void mult_simple_ptu22t(const MatrixA& a, const MatrixB& b, MatrixC& c)
{
    typedef typename MatrixC::value_type  value_type;
    const value_type z= math::zero(c[0][0]);    // if this are matrices we need their size
    
    set_to_zero(c);b

    // Temporary solution; dense matrices need to return const referencens
    MatrixA& aref= const_cast<MatrixA&>(a);
    MatrixB& bref= const_cast<MatrixB&>(b);

    size_t ari= &aref(1, 0) - &aref(0, 0), // how much is the offset of A's entry increased by incrementing row
	aci= &aref(0, 1) - &aref(0, 0), bri= &bref(1, 0) - &bref(0, 0), bci= &bref(0, 1) - &bref(0, 0);

    for (unsigned i= 0; i < c.num_rows(); i+=2)
	for (unsigned k= 0; k < c.num_cols(); k+=2) {
	    int ld= b.num_rows();
	    value_type tmp00= z, tmp01= z, tmp10= z, tmp11= z;

	    const value_type *begin_a= &aref[i][0], *end_a= &aref[i][a.num_cols()];
	    const value_type *begin_b= &bref[0][k];
	    for (; begin_a != end_a; begin_a+= aci, begin_b+= bri) {
		tmp00+= *begin_a * *begin_b;
		tmp01+= *begin_a * *(begin_b+bci);
		tmp10+= *(begin_a+ari) * *begin_b;
		tmp11+= *(begin_a+ari) * *(begin_b+bci);
	    }
	    assign::assign_sum::update(c[i][k], tmp00);
	    assign::assign_sum::update(c[i][k+1], tmp01);
	    assign::assign_sum::update(c[i+1][k], tmp10);
	    assign::assign_sum::update(c[i+1][k+1], tmp11);

#if 0
	    c[i][k]= tmp00; c[i][k+1]= tmp01;
	    c[i+1][k]= tmp10; c[i+1][k+1]= tmp11;
#endif
	}
}


// C must be square!
// no use of temporaries, let's see how the compiler handles this
void mult_simple_ptu22n(dense2D<double>& a, cm_type& b, dense2D<double>& c)
{
    for (unsigned i= 0; i < c.num_rows(); i+=2)
	for (unsigned k= 0; k < c.num_cols(); k+=2) {
	    int ld= b.num_rows();
	    double &tmp00= c[i][k], &tmp01= c[i][k+1], &tmp10=  c[i+1][k], &tmp11= c[i+1][k+1];
	    tmp00= 0.0, tmp01= 0.0, tmp10= 0.0, tmp11= 0.0;

	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()];
	    double *begin_b= &b[0][k];
	    for (; begin_a != end_a; ++begin_a, ++begin_b) {
		tmp00+= *begin_a * *begin_b;
		tmp01+= *begin_a * *(begin_b+ld);
		tmp10+= *(begin_a+ld) * *begin_b;
		tmp11+= *(begin_a+ld) * *(begin_b+ld);
	    }
	    //c[i][k]= tmp00; c[i][k+1]= tmp01;
	    //c[i+1][k]= tmp10; c[i+1][k+1]= tmp11;
	}
}


void mult_simple_ptu214(dense2D<double>& a, cm_type& b, dense2D<double>& c)
{

    // const double z= zero(
    for (unsigned i= 0; i < c.num_rows(); i+=2)
	for (unsigned k= 0; k < c.num_cols(); k+=4) {
	    int ld1= b.num_rows(), ld2= 2*ld1, ld3=3*ld1;
	    double tmp00= 0.0, tmp01= 0.0, tmp02= 0.0, tmp03= 0.0,
	  	   tmp10= 0.0, tmp11= 0.0, tmp12= 0.0, tmp13= 0.0;

	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()];
	    double *begin_b= &b[0][k];
	    for (; begin_a != end_a; ++begin_a, ++begin_b) {
		tmp00+= *begin_a * *begin_b;
		tmp01+= *begin_a * *(begin_b+ld1);
		tmp02+= *begin_a * *(begin_b+ld2);
		tmp03+= *begin_a * *(begin_b+ld3);
		tmp10+= *(begin_a+ld1) * *begin_b;
		tmp11+= *(begin_a+ld1) * *(begin_b+ld1);
		tmp12+= *(begin_a+ld1) * *(begin_b+ld2);
		tmp13+= *(begin_a+ld1) * *(begin_b+ld3);
	    }
	    c[i][k]= tmp00; c[i][k+1]= tmp01;
	    c[i][k+2]= tmp02; c[i][k+3]= tmp03;
	    c[i+1][k]= tmp10; c[i+1][k+1]= tmp11;
	    c[i+1][k+2]= tmp12; c[i+1][k+3]= tmp13;
	}
}


void mult_simple_ptu224(dense2D<double>& a, cm_type& b, dense2D<double>& c)
{
    for (unsigned i= 0; i < c.num_rows(); i+=2)
	for (unsigned k= 0; k < c.num_cols(); k+=4) {
	    int ld1= b.num_rows(), ld2= 2*ld1, ld3=3*ld1, lda1= a.num_rows();
	    double tmp000= 0.0, tmp001= 0.0, tmp002= 0.0, tmp003= 0.0,
	           tmp010= 0.0, tmp011= 0.0, tmp012= 0.0, tmp013= 0.0,
	  	   tmp100= 0.0, tmp101= 0.0, tmp102= 0.0, tmp103= 0.0,
	  	   tmp110= 0.0, tmp111= 0.0, tmp112= 0.0, tmp113= 0.0;

	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()];
	    double *begin_b= &b[0][k];
	    for (; begin_a != end_a; begin_a+= 2, begin_b+= 2) {
		tmp000+= *(begin_a) * *(begin_b);
		tmp001+= *(begin_a) * *(begin_b+ld1);
		tmp002+= *(begin_a) * *(begin_b+ld2);
		tmp003+= *(begin_a) * *(begin_b+ld3);
		tmp010+= *(begin_a+1) * *(begin_b+1);
		tmp011+= *(begin_a+1) * *(begin_b+1+ld1);
		tmp012+= *(begin_a+1) * *(begin_b+1+ld2);
		tmp013+= *(begin_a+1) * *(begin_b+1+ld3);
		tmp100+= *(begin_a+lda1) * *(begin_b);
		tmp101+= *(begin_a+lda1) * *(begin_b+ld1);
		tmp102+= *(begin_a+lda1) * *(begin_b+ld2);
		tmp103+= *(begin_a+lda1) * *(begin_b+ld3);
		tmp110+= *(begin_a+1+lda1) * *(begin_b+1);
		tmp111+= *(begin_a+1+lda1) * *(begin_b+1+ld1);
		tmp112+= *(begin_a+1+lda1) * *(begin_b+1+ld2);
		tmp113+= *(begin_a+1+lda1) * *(begin_b+1+ld3);
	    }
	    c[i][k]= tmp000 + tmp010; c[i][k+1]= tmp001 + tmp011;
	    c[i][k+2]= tmp002 + tmp012; c[i][k+3]= tmp003 + tmp013;
	    c[i+1][k]= tmp100 + tmp110; c[i+1][k+1]= tmp101 + tmp111;
	    c[i+1][k+2]= tmp102 + tmp112; c[i+1][k+3]= tmp103 + tmp113;
	}
}


void mult_simple_ptu224g(dense2D<double>& a, cm_type& b, dense2D<double>& c)
{
    using std::size_t;
    for (unsigned i= 0; i < c.num_rows(); i+=2)
	for (unsigned k= 0; k < c.num_cols(); k+=4) {
	    size_t ari= a.c_offset(1, 0), // how much is the offset of A's entry increased by incrementing row
		   aci= a.c_offset(0, 1), bri= b.c_offset(1, 0), bci= b.c_offset(0, 1);
	    double tmp000= 0.0, tmp001= 0.0, tmp002= 0.0, tmp003= 0.0,
	           tmp010= 0.0, tmp011= 0.0, tmp012= 0.0, tmp013= 0.0,
	  	   tmp100= 0.0, tmp101= 0.0, tmp102= 0.0, tmp103= 0.0,
	  	   tmp110= 0.0, tmp111= 0.0, tmp112= 0.0, tmp113= 0.0;

	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()];
	    double *begin_b= &b[0][k];
	    for (; begin_a != end_a; begin_a+= 2*aci, begin_b+= 2*bri) {
		tmp000+= *(begin_a+0*ari+0*aci) * *(begin_b+0*bri+0*bci);
		tmp001+= *(begin_a+0*ari+0*aci) * *(begin_b+0*bri+1*bci);
		tmp002+= *(begin_a+0*ari+0*aci) * *(begin_b+0*bri+2*bci);
		tmp003+= *(begin_a+0*ari+0*aci) * *(begin_b+0*bri+3*bci);
		tmp010+= *(begin_a+0*ari+1*aci) * *(begin_b+1*bri+0*bci);
		tmp011+= *(begin_a+0*ari+1*aci) * *(begin_b+1*bri+1*bci);
		tmp012+= *(begin_a+0*ari+1*aci) * *(begin_b+1*bri+2*bci);
		tmp013+= *(begin_a+0*ari+1*aci) * *(begin_b+1*bri+3*bci);
		tmp100+= *(begin_a+1*ari+0*aci) * *(begin_b+0*bri+0*bci);
		tmp101+= *(begin_a+1*ari+0*aci) * *(begin_b+0*bri+1*bci);
		tmp102+= *(begin_a+1*ari+0*aci) * *(begin_b+0*bri+2*bci);
		tmp103+= *(begin_a+1*ari+0*aci) * *(begin_b+0*bri+3*bci);
		tmp110+= *(begin_a+1*ari+1*aci) * *(begin_b+1*bri+0*bci);
		tmp111+= *(begin_a+1*ari+1*aci) * *(begin_b+1*bri+1*bci);
		tmp112+= *(begin_a+1*ari+1*aci) * *(begin_b+1*bri+2*bci);
		tmp113+= *(begin_a+1*ari+1*aci) * *(begin_b+1*bri+3*bci);
	    }
	    c[i][k]= tmp000 + tmp010; c[i][k+1]= tmp001 + tmp011;
	    c[i][k+2]= tmp002 + tmp012; c[i][k+3]= tmp003 + tmp013;
	    c[i+1][k]= tmp100 + tmp110; c[i+1][k+1]= tmp101 + tmp111;
	    c[i+1][k+2]= tmp102 + tmp112; c[i+1][k+3]= tmp103 + tmp113;
	}
}

// i, k, j == Outer, Middle, Inner
// Outer == row of A and C
// Middle == column of B and C
// Inner == column of A and row of B 
template <unsigned Outer, unsigned OuterMax, unsigned Middle, unsigned MiddleMax,
	  unsigned Inner, unsigned InnerMax>
struct double_matmat_mult_block 
    : public triple_unroll<Outer, OuterMax, Middle, MiddleMax, Inner, InnerMax>
{
    typedef triple_unroll<Outer, OuterMax, Middle, MiddleMax, Inner, InnerMax> base;
    typedef double_matmat_mult_block<base::next_outer, OuterMax, base::next_middle, MiddleMax, base::next_inner, InnerMax> next_t;
    typedef double       v_t;
    typedef std::size_t  s_t;

    void operator() (v_t &tmp00, v_t &tmp01, v_t &tmp02, v_t &tmp03, v_t &tmp04, 
		     v_t &tmp05, v_t &tmp06, v_t &tmp07, v_t &tmp08, v_t &tmp09, 
		     v_t &tmp10, v_t &tmp11, v_t &tmp12, v_t &tmp13, v_t &tmp14, v_t &tmp15, 
		     v_t *&begin_a, s_t &ari, s_t &aci, v_t *&begin_b, s_t &bri, s_t &bci)
    {
	tmp00+= begin_a[ this->outer * ari + this->inner * aci ] * begin_b[ this->inner * bri + this->middle * bci ];
	next_t()(tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
		 tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp00, 
		 begin_a, ari, aci, begin_b, bri, bci); 
    }

    template<typename MatrixC>
    void update(v_t &tmp00, v_t &tmp01, v_t &tmp02, v_t &tmp03, v_t &tmp04, 
		v_t &tmp05, v_t &tmp06, v_t &tmp07, v_t &tmp08, v_t &tmp09, 
		v_t &tmp10, v_t &tmp11, v_t &tmp12, v_t &tmp13, v_t &tmp14, v_t &tmp15, 
		MatrixC& c, s_t i, s_t k)
    {
	c[i + this->outer][k + this->middle]+= tmp00;
	next_t().update(tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
			tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp00, 
			c, i, k);
    }
};

template <unsigned OuterMax, unsigned MiddleMax, unsigned InnerMax>
struct double_matmat_mult_block<OuterMax, OuterMax, MiddleMax, MiddleMax, InnerMax, InnerMax>
    : public triple_unroll<OuterMax, OuterMax, MiddleMax, MiddleMax, InnerMax, InnerMax>
{
    typedef double       v_t;
    typedef std::size_t  s_t;

    void operator() (v_t &tmp00, v_t &, v_t &, v_t &, v_t &, 
		     v_t &, v_t &, v_t &, v_t &, v_t &, 
		     v_t &, v_t &, v_t &, v_t &, v_t &, v_t &, 
		     v_t *&begin_a, s_t &ari, s_t &aci, v_t *&begin_b, s_t &bri, s_t &bci)
    {
	tmp00+= begin_a[ this->outer * ari + this->inner * aci ] * begin_b[ this->inner * bri + this->middle * bci ];
    }

    template<typename MatrixC>
    void update(v_t &tmp00, v_t &, v_t &, v_t &, v_t &, 
		v_t &, v_t &, v_t &, v_t &, v_t &, 
		v_t &, v_t &, v_t &, v_t &, v_t &, v_t &, 
		MatrixC& c, s_t i, s_t k)
    {
	c[i + this->outer][k + this->middle]+= tmp00;
   }

};

template<unsigned Outer, unsigned Middle, unsigned Inner, typename MatrixA, typename MatrixB, typename MatrixC>
void double_matmat_mult_template(MatrixA& a, MatrixB& b, MatrixC& c)
{
    BOOST_STATIC_ASSERT(Outer * Middle * Inner <= 16);
 
    set_to_zero(c);
    double_matmat_mult_block<1, Outer, 1, Middle, 1, Inner> block;
    using std::size_t;
    for (size_t i= 0; i < c.num_rows(); i+= Outer)
	for (size_t k= 0; k < c.num_cols(); k+= Middle) {
	    size_t ari= a.c_offset(1, 0), // how much is the offset of A's entry increased by incrementing row
		   aci= a.c_offset(0, 1), bri= b.c_offset(1, 0), bci= b.c_offset(0, 1);
	    double tmp00= 0.0, tmp01= 0.0, tmp02= 0.0, tmp03= 0.0, tmp04= 0.0,
                   tmp05= 0.0, tmp06= 0.0, tmp07= 0.0, tmp08= 0.0, tmp09= 0.0,
 		   tmp10= 0.0, tmp11= 0.0, tmp12= 0.0, tmp13= 0.0, tmp14= 0.0, tmp15= 0.0;
	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()], *begin_b= &b[0][k];

	    for (; begin_a != end_a; begin_a+= Inner*aci, begin_b+= Inner*bri)
		block(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
		      tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, 
		      begin_a, ari, aci, begin_b, bri, bci); 
	    block.update(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
			 tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, 
			 c, i, k);
	}
}


// i, k, j == Outer, Inner, __
// i == Outer == row of A and C
// k == Inner == column of B and C
// j == column of A and row of B  is not unrolled !!!
template <unsigned Outer, unsigned OuterMax, unsigned Inner, unsigned InnerMax>
struct twice_double_matmat_mult_block 
    : public double_unroll<Outer, OuterMax, Inner, InnerMax>
{
    typedef double_unroll<Outer, OuterMax, Inner, InnerMax> base;
    typedef twice_double_matmat_mult_block<base::next_outer, OuterMax, base::next_inner, InnerMax> next_t;
    typedef double       v_t;
    typedef std::size_t  s_t;

    void operator() (v_t &tmp00, v_t &tmp01, v_t &tmp02, v_t &tmp03, v_t &tmp04, 
		     v_t &tmp05, v_t &tmp06, v_t &tmp07, v_t &tmp08, v_t &tmp09, 
		     v_t &tmp10, v_t &tmp11, v_t &tmp12, v_t &tmp13, v_t &tmp14, v_t &tmp15, 
		     v_t *begin_a, s_t &ari, v_t *begin_b, s_t &bci)
    {
	tmp00+= begin_a[ this->outer * ari ] * begin_b[ this->inner * bci ];
	next_t()(tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
		 tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp00, 
		 begin_a, ari, begin_b, bci); 
    }

    template<typename MatrixC>
    void update(v_t &tmp00, v_t &tmp01, v_t &tmp02, v_t &tmp03, v_t &tmp04, 
		v_t &tmp05, v_t &tmp06, v_t &tmp07, v_t &tmp08, v_t &tmp09, 
		v_t &tmp10, v_t &tmp11, v_t &tmp12, v_t &tmp13, v_t &tmp14, v_t &tmp15, 
		MatrixC& c, s_t i, s_t k)
    {
	c[i + this->outer][k + this->inner]+= tmp00;
	next_t().update(tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
			tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp00, 
			c, i, k);
    }
};

template <unsigned OuterMax, unsigned InnerMax>
struct twice_double_matmat_mult_block<OuterMax, OuterMax, InnerMax, InnerMax>
    : public double_unroll<OuterMax, OuterMax, InnerMax, InnerMax>
{
    typedef double       v_t;
    typedef std::size_t  s_t;

    void operator() (v_t &tmp00, v_t &, v_t &, v_t &, v_t &, 
		     v_t &, v_t &, v_t &, v_t &, v_t &, 
		     v_t &, v_t &, v_t &, v_t &, v_t &, v_t &, 
		     v_t *begin_a, s_t &ari, v_t *begin_b, s_t &bci)
    {
	tmp00+= begin_a[ this->outer * ari ] * begin_b[ this->inner * bci ];
    }

    template<typename MatrixC>
    void update(v_t &tmp00, v_t &, v_t &, v_t &, v_t &, 
		v_t &, v_t &, v_t &, v_t &, v_t &, 
		v_t &, v_t &, v_t &, v_t &, v_t &, v_t &, 
		MatrixC& c, s_t i, s_t k)
    {
	c[i + this->outer][k + this->inner]+= tmp00;
   }

};


template<unsigned Outer, unsigned Inner, typename MatrixA, typename MatrixB, typename MatrixC>
void twice_double_matmat_mult_template(MatrixA& a, MatrixB& b, MatrixC& c)
{
    BOOST_STATIC_ASSERT(Outer * Inner <= 16);
 
    set_to_zero(c);
    twice_double_matmat_mult_block<1, Outer, 1, Inner> block;
    using std::size_t;
    for (size_t i= 0; i < c.num_rows(); i+= Outer)
	for (size_t k= 0; k < c.num_cols(); k+= Inner) {
	    size_t ari= a.c_offset(1, 0), // how much is the offset of A's entry increased by incrementing row
		   aci= a.c_offset(0, 1), bri= b.c_offset(1, 0), bci= b.c_offset(0, 1);
	    double tmp00= 0.0, tmp01= 0.0, tmp02= 0.0, tmp03= 0.0, tmp04= 0.0,
                   tmp05= 0.0, tmp06= 0.0, tmp07= 0.0, tmp08= 0.0, tmp09= 0.0,
 		   tmp10= 0.0, tmp11= 0.0, tmp12= 0.0, tmp13= 0.0, tmp14= 0.0, tmp15= 0.0;
	    double *begin_a= &a[i][0], *end_a= &a[i][a.num_cols()], *begin_b= &b[0][k];

	    for (; begin_a != end_a; begin_a+= aci, begin_b+= bri)
		block(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
		      tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, 
		      begin_a, ari, begin_b, bci); 
	    block.update(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
			 tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, 
			 c, i, k);
	}
}



int test_main(int argc, char* argv[])
{
    //double_unroll<1, 4, 1, 6>()();
    //triple_wrapper<1, 4, 1, 4, 1, 4>() ();
    //triple_unroll<4, 4, 4, 4, 4, 4>() ();
    //triple_unroll<3, 3, 6, 6, 4, 4>() ();

    unsigned steps= 32, max_size= 128, size= 32; 
    if (argc > 2) {
	steps= atoi(argv[1]); max_size= atoi(argv[2]);
    }
    
    mtl::dense2D<double>               da(size, size), db(size, size), dc(size, size);
    cm_type                            dat(size, size), dbt(size, size), dct(size, size);
    hessian_setup(da, 1.0);
    hessian_setup(db, 2.0); 
    hessian_setup(dbt, 2.0); 

    time_series(da, dbt, dc, mult_simple_ptu22t<rm_type, cm_type, rm_type>, "Same templated", steps, max_size);
    time_series(da, dbt, dc, mult_simple_ptu22, "Simple mult (pointers trans unrolled 2x1x2)", steps, max_size);
    return 0;

    time_series(da, dbt, dc, twice_double_matmat_mult_template<1, 1, rm_type, cm_type, rm_type>, 
		"Simple mult (pointers trans unrolled 1x1 template)", steps, max_size);
    time_series(da, dbt, dc, twice_double_matmat_mult_template<2, 2, rm_type, cm_type, rm_type>, 
		"Simple mult (pointers trans unrolled 2x2 template)", steps, max_size);
    time_series(da, dbt, dc, twice_double_matmat_mult_template<2, 4, rm_type, cm_type, rm_type>, 
		"Simple mult (pointers trans unrolled 2x4 template)", steps, max_size);
    time_series(da, dbt, dc, twice_double_matmat_mult_template<4, 2, rm_type, cm_type, rm_type>, 
		"Simple mult (pointers trans unrolled 4x2 template)", steps, max_size);
    time_series(da, dbt, dc, twice_double_matmat_mult_template<4, 4, rm_type, cm_type, rm_type>, 
		"Simple mult (pointers trans unrolled 4x4 template)", steps, max_size);


    time_series(da, dbt, dc, double_matmat_mult_template<2, 2, 2, rm_type, cm_type, rm_type>, 
		"Simple mult (pointers trans unrolled 2x2x2 template)", steps, max_size);
    time_series(da, dbt, dc, double_matmat_mult_template<1, 1, 1, rm_type, cm_type, rm_type>, 
		"Simple mult (pointers trans unrolled 1x1x1 template)", steps, max_size);
    time_series(da, dbt, dc, double_matmat_mult_template<2, 4, 1, rm_type, cm_type, rm_type>, 
		"Simple mult (pointers trans unrolled 2x1x4 template)", steps, max_size);
    time_series(da, db, dc, double_matmat_mult_template<2, 4, 1, rm_type, rm_type, rm_type>, 
		"Simple mult (pointers unrolled 2x1x4 template)", steps, max_size);
    time_series(da, dbt, dc, double_matmat_mult_template<2, 4, 2, rm_type, cm_type, rm_type>, 
		"Simple mult (pointers trans unrolled 2x2x4 template)", steps, max_size);

    time_series(da, dbt, dc, mult_simple_ptu224g, "Simple mult (pointers trans unrolled 2x2x4 generic)", steps, max_size);
    time_series(da, dbt, dc, mult_simple_ptu224, "Simple mult (pointers trans unrolled 2x2x4)", steps, max_size);
    time_series(da, dbt, dc, mult_simple_ptu214, "Simple mult (pointers trans unrolled 2x1x4)", steps, max_size);
    time_series(da, dbt, dc, mult_simple_ptu22n, "Simple mult (pointers trans unrolled 2x1x2 no temps)", steps, max_size);
    time_series(da, dbt, dc, mult_simple_ptu22, "Simple mult (pointers trans unrolled 2x1x2)", steps, max_size);
    time_series(da, dbt, dc, mult_simple_ptu4, "Simple mult (pointers trans unrolled 4)", steps, max_size);
    time_series(da, dbt, dc, mult_simple_ptu, "Simple mult (pointers trans unrolled 2)", steps, max_size);
    time_series(da, dbt, dc, mult_simple_pt, "Simple mult (pointers transposed)", steps, max_size);
    time_series(da, db, dc, mult_simple_p, "Simple mult (pointers)", steps, max_size);
    time_series(da, db, dc, mult_simple, "Simple mult", steps, max_size);

    return 0;
}
