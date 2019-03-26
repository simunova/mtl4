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
#include <algorithm>
#include <cmath>
#include <string>

#include <boost/numeric/mtl/mtl.hpp> 

using namespace std;  


typedef mtl::dense2D<int> matrix_type;
using mtl::irange; using mtl::iall;


matrix_type inline my_row(int i, int, matrix_type& A)
{
    return A[irange(i, i+1)][iall];
}

matrix_type inline my_col(int, int j, matrix_type& A)
{
    return A[iall][irange(j, j+1)];
}

matrix_type inline my_block(int i, int j, matrix_type& A)
{
    int ib= i/3*3, jb= j/3*3;
    return A[irange(ib, ib+3)][irange(jb, jb+3)];
}

int inline nzero(const matrix_type& A)
{
    int n= 0;
    for (unsigned i= 0; i < num_rows(A); i++)
	for (unsigned j= 0; j < num_cols(A); j++)
	    if (A[i][j]) n++;
    return n;
}

// My relevant sub-matrices
struct my_sub_t
{
    my_sub_t(int i, int j, matrix_type& A) 
	: i(i), j(j), mr(my_row(i, j, A)), mc(my_col(i, j, A)), mb(my_block(i, j, A))
    {}

    int         i, j;
    matrix_type mr, mc, mb;
};

bool inline conflict(int v, const matrix_type& A)
{
    for (unsigned i= 0; i < num_rows(A); i++)
	for (unsigned j= 0; j < num_cols(A); j++)
	    if (A[i][j] == v) return true;
    return false;
}

bool inline conflict(int v, const my_sub_t& sub)
{
    return conflict(v, sub.mr) || conflict(v, sub.mc) || conflict(v, sub.mb);
}


struct entry
{
    entry(int i, int j, matrix_type& A) 
	: i(i), j(j), nnz(nzero(my_row(i, j, A)) + nzero(my_col(i, j, A)) + nzero(my_block(i, j, A)))
    {}

    friend inline std::ostream& operator<< (std::ostream& stream, const entry& e) 
    {
	return stream << "[" << e.i << ", " << e.j << " = " << e.nnz << "]";
    }

    bool operator<(const entry& other) const { return nnz > other.nnz; } // to sort 

    int         i, j, nnz;
};

template <typename T>
inline std::ostream& operator<< (std::ostream& stream, const std::vector<T>& v) 
{
    stream << "(";
    for (typename std::vector<T>::const_iterator it= v.begin(); it != v.end(); ++it)
	stream << *it << ", ";
    return stream << ")";
}


// void inline wait() { char c; std::cin >> c; }


void solve(const char* file_name)
{
    mtl::dense2D<int> A;
    mtl::io::matrix_market_istream(file_name) >> A;
    cout << "Read from " << file_name <<  " is \n"  << A;

    std::vector<entry> to_fill;
    for (int i= 0; i < 9; i++)
	for (int j= 0; j < 9; j++)
	    if (A[i][j] == 0)
		to_fill.push_back(entry(i, j, A));
    sort(to_fill.begin(), to_fill.end());

    for (unsigned pos= 0;;) {
	int i= to_fill[pos].i, j= to_fill[pos].j;
	my_sub_t my_sub(i, j, A);
	int v= A[i][j] + 1;
	while (v < 10 && conflict(v, my_sub))
	    v++;

	// cout << "Matrix is now: \n" << A;
	if (v == 10) {// nothing works -> back track
	    if (pos-- == 0) 
		throw "No solution found";
	    A[i][j]= 0;
	} else {
	    A[i][j]= v;
	    if (++pos == to_fill.size()) {
		cout << "Solution is:\n" << A;
		return;
	    }
	}
    }

}


int main(int argc, char* argv[])
{
    solve(argc > 1 ? argv[1] : "matrix_market/sudoku_easy.mtx");
    return 0;
}
