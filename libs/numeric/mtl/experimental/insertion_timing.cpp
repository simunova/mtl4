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
// #include <boost/test/minimal.hpp>

#include <boost/numeric/mtl/operation/update.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>

using namespace mtl;
using namespace std;
using namespace boost;

template <typename Generator>
double generator_overhead(int s, Generator gen)
{    
    boost::timer t;
    while (! gen.finished() ) {
	int      r, c;
	double   v;
	gen(r, c, v);
    }
    return t.elapsed();
}

template <typename Generator>
double insert_petsc(int s, Generator gen, double overhead)
{
    return 0.0;
}

template <typename Generator>
double insert_mtl2(int s, Generator gen, double overhead)
{
    return 0.0;
}

template <typename Generator>
double insert_mtl4(int s, Generator gen, double overhead)
{
    typedef typename Generator::value_type                                             value_type;
    typedef mat::parameters<row_major, mtl::index::c_index, non_fixed::dimensions>   parameters;
    typedef compressed2D<value_type, parameters>                                       matrix_type;
    matrix_type   matrix(non_fixed::dimensions(s, s)); 
  
    boost::timer t;
    compressed2D_inserter<value_type, parameters>  inserter(matrix, 7);
    while (! gen.finished() ) {
	int      r, c;
	double   v;
	gen(r, c, v);
	inserter(r, c) << v;
	// cout << "A[" << r << ", " << c << "] = " << v << '\n';
    }
    return t.elapsed();
}

// produce a matrix with s nonzeros, randomly distributed
template <typename Size, typename Value, typename RandomGen,
	  typename RowDistribution = uniform_int<Size>,
	  typename ColumnDistribution = RowDistribution>
struct random_generator
{
    typedef Value value_type;

    explicit random_generator( Size s, RandomGen mygen = RandomGen() ) :
	mygen( mygen ), nnz( 5 * s ), 
	pick_row( RowDistribution(0,s-1) ), pick_col( ColumnDistribution(0,s-1) ) {}
	
    bool finished() 
    {
	return !(bool)nnz;
    }
    
    void operator()( Size& i, Size& j, value_type& v )
    {
	i = pick_row( mygen );
	j = pick_col( mygen );
	v = 1.;
	nnz--;
    }

    RandomGen          mygen;
    Size               nnz;
    RowDistribution    pick_row; 
    ColumnDistribution pick_col;
};
    
// produce a matrix with s nonzeros, randomly distributed
template <typename Size, typename Value, typename RandomGen,
	  typename ColumnDistribution = uniform_int<Size> >
struct random_column_generator
{
    typedef Value value_type;

    explicit random_column_generator( Size s, RandomGen mygen = RandomGen() ) :
	mygen( mygen ), s( s ), nnz( 5 * s ), 
	pick_col( ColumnDistribution(0,s-1) ) {}
	
    bool finished() 
    {
	return !(bool)nnz;
    }
    
    void operator()( Size& i, Size& j, value_type& v )
    {
	i = s - (nnz+4) / 5;
	// if (i < 0 || i >= s) cout << "i = " << i << ", s = " << s << ", nnz = " << nnz << '\n';
	assert(i >= 0 && i < s);
	j = pick_col( mygen );
	assert(j >= 0 && j < s);
	v = 1.;
	nnz--; 
    }

    RandomGen          mygen;
    Size               s, nnz;
    ColumnDistribution pick_col;
};
    

template <typename Size, typename Value>
struct poisson_generator
{
    typedef Value  value_type;

    explicit poisson_generator(int s) : s(s), count(0), row(0), offset(2) {     // s must be 2^k 100
        if (s == 9) {                                      // only to test a 9x9 matrix
	    d1= d2= 3;                                     // only to test a 9x9 matrix
	} else {
	    assert(s % 100 == 0);
	    d1= 100, d2= s/d1;
	    for (; d2 > d1; d1<<= 1) d2= s / d1;
	}
	nnz= 5 * s - 2 * d1 - 2 * d2;
    }
    
    bool finished() { 
	return count >= nnz; 
    }

    bool valid_offset() {
	switch (offset) {
	    case 0: return row >= d2; // northern nb
	    case 1: return row % d2;  // western nb
	    case 2: return true;
	    case 3: return (row + 1) % d2; // eastern nb
	    case 4: return row + d2 < s;   // southern nb
	}
	assert(true); // shouldn't be reached
	return false;
    }

    void next_offset() {
	offset++;
	if (offset == 5) offset= 0, row++;
    }

    void next_nnz() {
	do 
	    {next_offset(); } 
	while (! valid_offset());
	count++;
    }

    void operator() (Size& r, Size& c, Value& v) {
	r= row;
	switch (offset) {
	    case 0: c= r - d2; v= -1; break;
	    case 1: c= r - 1; v= -1; break;
	    case 2: c= r; v= 4; break;
	    case 3: c= r + 1; v= -1; break;
	    case 4: c= r + d2; v= -1; 
	}
	next_nnz();
    }	    

    int s, d1, d2, count, 
	nnz,                   // number of non-zeros in matrix
	row, offset;           // current row and which entry in the row
};



template <typename Generator>
void run(int max_size)
{
    for (int s= 10000; s < max_size; s*= 2) {
	cout << "Size is " << s << '\n';

	Generator gen0(s);
	double overhead= generator_overhead(s, gen0);
	cout << "Generator overhead: " << overhead << "s\n";

	Generator gen1(s);
	insert_petsc(s, gen1, overhead);

	Generator gen2(s);
	insert_mtl2(s, gen2, overhead);
	
	Generator gen3(s);
	double mtl4t= insert_mtl4(s, gen3, overhead);
	cout << "Inserting into MTL 4: " << mtl4t << "s\n";
    }
}

void check_dims()
{
    for (int s= 100;  s < 1000000; s*= 2) {
	int d1= 10, d2= s/d1;
	for (; d2 > d1; d1<<= 1) 
	    d2= s / d1;
	std::cout << d1 << " * " << d2 << " = " << s << '\n';
    }
}

int main(int argc, char* argv[])
{
    // check_dims();
    poisson_generator<int, double> poisson_9(9);
    insert_mtl4(9, poisson_9, 0.0);

    minstd_rand gen;
    random_generator<int, double, minstd_rand> random_9(9);
    insert_mtl4(9, random_9, 0.0);

    if( argc > 1 ) {
	run<poisson_generator<int, double> > (atoi(argv[1]));
	run<random_generator<int, double, minstd_rand> > (atoi(argv[1]));
	run<random_column_generator<int, double, minstd_rand> > (atoi(argv[1]));
    }

    return 0;
}
