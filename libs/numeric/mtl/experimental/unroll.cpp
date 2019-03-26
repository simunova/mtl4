// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschr채nkt), www.simunova.com.
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

using namespace std;

int const vector_size = 1000, // 1000
          repetitions = 500000;
 
vector<double> gv1(vector_size, 2.0), gv2(vector_size, 3.0),
               gv3(vector_size, 4.0), gv4(vector_size, 5.0),
               gv5(vector_size, 6.0), gv6(vector_size, 7.0);

template <typename F>
void time_dot(std::string fname, F f)
{

    boost::timer start;
    double result;
    for (int i= 0; i < repetitions; i++) // 1000000
	result= f(gv1, gv2);
    double duration = start.elapsed();
    cout << fname << ": " << duration / repetitions * 1000000 << "탎" 
      // << ", result = " << result 
	 << "\n";
}

inline double 
dot(vector<double> const& v1, vector<double> const& v2)
{
    double sum= 0.0;
    for (unsigned i= 0; i < v1.size(); i++)
	sum+= v1[i] * v2[i];

    return sum;
}

inline double 
dot2(vector<double> const& v1, vector<double> const& v2)
{
    double sum= 0.0, sum2 = 0.0;
    for (unsigned i= 0; i < v1.size(); i+= 2) {
	sum+= v1[i] * v2[i];
	sum2+= v1[i+1] * v2[i+1];
    }

    return sum + sum2;
}

struct dot2slow_t
{
    double sum;
    struct {
	double sum2;
    } xx;
};

double dot2slow(vector<double> const& v1, vector<double> const& v2)
{
    dot2slow_t data;

    data.sum= 0.0, data.xx.sum2 = 0.0;
    for (unsigned i= 0; i < v1.size(); i+= 2) {
	data.sum+= v1[i] * v2[i];
	data.xx.sum2+= v1[i+1] * v2[i+1];
    }

    return data.sum + data.xx.sum2;
}

double dot4(vector<double> const& v1, vector<double> const& v2)
{
    double sum= 0.0, sum2 = 0.0, sum3= 0.0, sum4= 0.0;
    for (unsigned i= 0; i < v1.size(); i+= 4) {
	sum+= v1[i] * v2[i];
	sum2+= v1[i+1] * v2[i+1];
	sum3+= v1[i+2] * v2[i+2];
	sum4+= v1[i+3] * v2[i+3];
    }

    return sum + sum2 + sum3 + sum4;
}

template <unsigned Depth, unsigned MaxDepth>
struct dot_block
{
    static unsigned const offset= MaxDepth - Depth;

    void operator() (vector<double> const& v1, vector<double> const& v2, unsigned i,
		     double& s0, double& s1, double& s2, double& s3,
		     double& s4, double& s5, double& s6, double& s7)
    {
	s0+= v1[ i + offset ] * v2[ i + offset ];
	dot_block<Depth-1, MaxDepth>() (v1, v2, i, s1, s2, s3, s4, s5, s6, s7, s0);
    }
};

template <unsigned MaxDepth>
struct dot_block<1, MaxDepth>
{
    static unsigned const offset= MaxDepth - 1;

    void operator() (vector<double> const& v1, vector<double> const& v2, unsigned i,
		     double& s0, double&, double&, double&,
		     double&, double&, double&, double&)
    {
	s0+= v1[ i + offset ] * v2[ i + offset ];
    }
};


template <unsigned Depth>
double unrolled_dot(vector<double> const& v1, vector<double> const& v2)
{
    // check v1.size() == v2.size();
    unsigned size= v1.size(), blocks= size / Depth, blocked_size= blocks * Depth;

    double s0= 0.0, s1= 0.0, s2= 0.0, s3= 0.0, s4= 0.0, s5= 0.0, s6= 0.0, s7= 0.0;
    for (unsigned i= 0; i < blocked_size; i+= Depth)
	dot_block<Depth, Depth>()(v1, v2, i, s0, s1, s2, s3, s4, s5, s6, s7);

    s0+= s1 + s2 + s3 + s4 + s5 + s6 + s7;
    for (unsigned i= blocked_size; i < size; ++i)
	s0+= v1[i] * v2[i];
    return s0;
}


template <unsigned Depth>
struct sum_helper_t
{
    double operator()(double& s0, double& s1, double& s2, double& s3,
		      double& s4, double& s5, double& s6, double& s7)
    {
	return s0 + sum_helper_t<Depth-1>()(s1, s2, s3, s4, s5, s6, s7, s0);
    }
};

template <>
struct sum_helper_t<1>
{
    double operator()(double& s0, double&, double&, double&,
		      double&, double&, double&, double&)
    {
	return s0;
    }
};



template <unsigned Depth>
struct unrolled_dot_helper_t
{
    double inline operator()(vector<double> const& v1, vector<double> const& v2,
		      double& s0, double& s1, double& s2, double& s3,
		      double& s4, double& s5, double& s6, double& s7)
    {
	unsigned size= v1.size(), blocks= size / Depth, blocked_size= blocks * Depth;
	for (unsigned i= 0; i < blocked_size; i+= Depth)
	    dot_block<Depth, Depth>()(v1, v2, i, s0, s1, s2, s3, s4, s5, s6, s7);
	for (unsigned i= blocked_size; i < size; ++i)
	    s0+= v1[i] * v2[i];
	return sum_helper_t<Depth>()(s0, s1, s2, s3, s4, s5, s6, s7);
    }
};


// #define TRICKY

template <unsigned Depth>
struct unrolled_dot_t
{
    double inline operator()(vector<double> const& v1, vector<double> const& v2)
    {
#ifdef TRICKY
	double s0= 0.0, s1= 0.0, s2= 0.0, s3= 0.0, s4= 0.0, s5= 0.0, s6= 0.0, s7= 0.0;
	return unrolled_dot_helper_t<Depth>()(v1, v2, s0, s1, s2, s3, s4, s5, s6, s7);
#else
	unsigned size= v1.size(), blocks= size / Depth, blocked_size= blocks * Depth;

	double s0= 0.0, s1= 0.0, s2= 0.0, s3= 0.0, s4= 0.0, s5= 0.0, s6= 0.0, s7= 0.0;
	for (unsigned i= 0; i < blocked_size; i+= Depth)
	    dot_block<Depth, Depth>()(v1, v2, i, s0, s1, s2, s3, s4, s5, s6, s7);
	
	s0+= s1 + s2 + s3 + s4 + s5 + s6 + s7;
	for (unsigned i= blocked_size; i < size; ++i)
	    s0+= v1[i] * v2[i];
	return s0;
#endif
    }
};
	
template <>
struct unrolled_dot_t<4>
{
    double inline operator()(vector<double> const& v1, vector<double> const& v2)
    {
#ifdef TRICKY
	double s0= 0.0, s1= 0.0, s2= 0.0, s3= 0.0;
	return unrolled_dot_helper_t<4>()(v1, v2, s0, s1, s2, s3, s0, s1, s2, s3);
#else
	static const unsigned Depth= 4;
	unsigned size= v1.size(), blocks= size / Depth, blocked_size= blocks * Depth;

	double s0= 0.0, s1= 0.0, s2= 0.0, s3= 0.0;
	for (unsigned i= 0; i < blocked_size; i+= Depth)
	    dot_block<Depth, Depth>()(v1, v2, i, s0, s1, s2, s3, s0, s1, s2, s3);
	
	s0+= s1 + s2 + s3;
	for (unsigned i= blocked_size; i < size; ++i)
	    s0+= v1[i] * v2[i];
	return s0;
#endif
    }
};

template <unsigned Depth>
double unrolled_dot_id(vector<double> const& v1, vector<double> const& v2)
{
    return unrolled_dot_t<Depth>()(v1, v2);
}



// ===========================
// Enough dot product !!!!!!!!
// ===========================
	


inline void
expr(vector<double>& v1, vector<double> const& v2,
     vector<double> const& v3, vector<double> const& v4,
     double s1, double s2, double s3)
{
    for (unsigned i= 0; i < v1.size(); i++)
	v1[i] = s1 * v2[i] + s2 * v3[i] + s3 * v4[i];
}


inline void
expr4(vector<double>& v1, vector<double> const& v2,
     vector<double> const& v3, vector<double> const& v4,
     double s1, double s2, double s3)
{
    double t0= 0.0, t1= 0.0, t2= 0.0, t3= 0.0;
    for (unsigned i= 0; i < v1.size(); i+= 4) {
	t0=  s1 * v2[i] + s2 * v3[i] + s3 * v4[i];
	t1=  s1 * v2[i+1] + s2 * v3[i+1] + s3 * v4[i+1];
	t2=  s1 * v2[i+2] + s2 * v3[i+2] + s3 * v4[i+2];
	t3=  s1 * v2[i+3] + s2 * v3[i+3] + s3 * v4[i+3];
	v1[i]= t0; v1[i+1]= t1; v1[i+2]= t2; v1[i+3]= t3;
    }
}


inline void
expr8(vector<double>& v1, vector<double> const& v2,
     vector<double> const& v3, vector<double> const& v4,
     double s1, double s2, double s3)
{
    double t0= 0.0, t1= 0.0, t2= 0.0, t3= 0.0, t4= 0.0, t5= 0.0, t6= 0.0, t7= 0.0;
    for (unsigned i= 0; i < v1.size(); i+= 8) {
	t0=  s1 * v2[i] + s2 * v3[i] + s3 * v4[i];
	t1=  s1 * v2[i+1] + s2 * v3[i+1] + s3 * v4[i+1];
	t2=  s1 * v2[i+2] + s2 * v3[i+2] + s3 * v4[i+2];
	t3=  s1 * v2[i+3] + s2 * v3[i+3] + s3 * v4[i+3];
	t4=  s1 * v2[i+4] + s2 * v3[i+4] + s3 * v4[i+4];
	t5=  s1 * v2[i+5] + s2 * v3[i+5] + s3 * v4[i+5];
	t6=  s1 * v2[i+6] + s2 * v3[i+6] + s3 * v4[i+6];
	t7=  s1 * v2[i+7] + s2 * v3[i+7] + s3 * v4[i+7];
	v1[i]= t0; v1[i+1]= t1; v1[i+2]= t2; v1[i+3]= t3;
	v1[i+4]= t4; v1[i+5]= t5; v1[i+6]= t6; v1[i+7]= t7;
    }
}


template <typename F>
void time_expr(std::string fname, F f)
{

    boost::timer start;
    for (int i= 0; i < repetitions; i++) // 1000000
	f(gv1, gv2, gv3, gv4, 2.0, 3.0, 4.0); 
    double duration = start.elapsed();
    cout << fname << ": " << duration / repetitions * 1000000 << "탎" 
      // << ", result = " << result 
	 << "\n";
}

// =============================
// Again with binary combination
// =============================





inline void
bexpr(vector<double>& v1, vector<double> const& v2,
     vector<double> const& v3, 
     double s1, double s2)
{
    for (unsigned i= 0; i < v1.size(); i++)
	v1[i] = s1 * v2[i] + s2 * v3[i];
}


inline void
bexpr4(vector<double>& v1, vector<double> const& v2,
     vector<double> const& v3, 
     double s1, double s2)
{
    double t0= 0.0, t1= 0.0, t2= 0.0, t3= 0.0;
    for (unsigned i= 0; i < v1.size(); i+= 4) {
	t0=  s1 * v2[i] + s2 * v3[i];
	t1=  s1 * v2[i+1] + s2 * v3[i+1];
	t2=  s1 * v2[i+2] + s2 * v3[i+2];
	t3=  s1 * v2[i+3] + s2 * v3[i+3];
	v1[i]= t0; v1[i+1]= t1; v1[i+2]= t2; v1[i+3]= t3;
    }
}


inline void
bexpr8(vector<double>& v1, vector<double> const& v2,
     vector<double> const& v3, 
     double s1, double s2)
{
    double t0= 0.0, t1= 0.0, t2= 0.0, t3= 0.0, t4= 0.0, t5= 0.0, t6= 0.0, t7= 0.0;
    for (unsigned i= 0; i < v1.size(); i+= 8) {
	t0=  s1 * v2[i] + s2 * v3[i];
	t1=  s1 * v2[i+1] + s2 * v3[i+1];
	t2=  s1 * v2[i+2] + s2 * v3[i+2];
	t3=  s1 * v2[i+3] + s2 * v3[i+3];
	t4=  s1 * v2[i+4] + s2 * v3[i+4];
	t5=  s1 * v2[i+5] + s2 * v3[i+5];
	t6=  s1 * v2[i+6] + s2 * v3[i+6];
	t7=  s1 * v2[i+7] + s2 * v3[i+7];
	v1[i]= t0; v1[i+1]= t1; v1[i+2]= t2; v1[i+3]= t3;
	v1[i+4]= t4; v1[i+5]= t5; v1[i+6]= t6; v1[i+7]= t7;
    }
}


template <typename F>
void time_bexpr(std::string fname, F f)
{

    boost::timer start;
    for (int i= 0; i < repetitions; i++) // 1000000
	f(gv1, gv2, gv3, 2.0, 3.0); 
    double duration = start.elapsed();
    cout << fname << ": " << duration / repetitions * 1000000 << "탎" 
      // << ", result = " << result 
	 << "\n";
}




int test_main(int argc, char* argv[])
{
    time_bexpr("init      ", bexpr);
    time_bexpr("reg bexpr  ", bexpr);
    time_bexpr("bexpr 4    ", bexpr4);
    time_bexpr("bexpr 8    ", bexpr8);

    time_expr("reg expr  ", expr);
    time_expr("expr 4    ", expr4);
    time_expr("expr 8    ", expr8);

    time_dot("init      ", dot);
    time_dot("regular   ", dot);
    time_dot("unrolled 2", dot2);
    time_dot("struct 2  ", dot2slow);
    time_dot("unrolled 4", dot4);

    cout << "--------------------\n";
    time_dot("template 1", unrolled_dot_id<1>);
    time_dot("template 2", unrolled_dot_id<2>);
    time_dot("template 3", unrolled_dot_id<3>);
    time_dot("template 4", unrolled_dot_id<4>);
    
    cout << "--------------------\n";
    time_dot("template 1", unrolled_dot<1>);
    time_dot("template 2", unrolled_dot<2>);
    time_dot("template 3", unrolled_dot<3>);
    time_dot("template 4", unrolled_dot<4>);
    time_dot("template 5", unrolled_dot<5>);
    time_dot("template 6", unrolled_dot<6>);
    time_dot("template 7", unrolled_dot<7>);
    time_dot("template 8", unrolled_dot<8>);
    time_dot("template 9", unrolled_dot<9>);
    time_dot("template10", unrolled_dot<10>);
    time_dot("template11", unrolled_dot<11>);
    time_dot("template12", unrolled_dot<12>);
    time_dot("template13", unrolled_dot<13>);
    time_dot("template14", unrolled_dot<14>);
    time_dot("template15", unrolled_dot<15>);
    time_dot("template16", unrolled_dot<16>);

    return 0;
}
