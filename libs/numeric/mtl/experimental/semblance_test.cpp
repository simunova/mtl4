// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschrÃ¤nkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;
using namespace mtl;

template <typename Matrix>
void init_data(Matrix& T)
{
    for (unsigned i= 0; i < num_rows(T); i++)
	for (unsigned j= 0; j < num_cols(T); j++)
	    T[i][j]= i + j + 2; // add 2 because i and j are zero-based
}

template <typename Matrix, typename Vector>
void w_semblance(const Matrix& T, Vector& semb, int nwin, int lwin, int linc)
{
    int nsamp= num_rows(T), ntr= num_cols(T);
    dense_vector<double> sumsqr(nsamp), sqrsum(nsamp);

    for (int i= 0; i < nsamp; i++) {
	sumsqr[i]= unary_dot(T[i][iall]);
	sqrsum[i]= square(sum(T[i][iall]));
    }
    
    for (int i= 1, ll= 0; i <= nwin; i++, ll+= linc) {
	int kkend= ll + kkend > nsamp ? nsamp - ll : lwin;
	irange r(ll, kkend);
	double sumsq= sum(sumsqr[r]),
	       sumampsq= sum(sqrsum[r]),
	       value= ntr * sumsq;
	semb[i]= value != 0.0 ? sumampsq / value : 0.0;
    }
}


int main() 
{
    //int nsamples= 30001, ntrc= 2000;
    int nsamples= 8, ntrc= 6;

    dense2D<float> traces(nsamples, ntrc);
    init_data(traces);
    // cout << "T is\n" << traces;

    dense_vector<double> semb(nsamples/2 + 1, 0.0);
    int lwin= 5, linc= 2, nwin= size(semb) - linc;
    w_semblance(traces, semb, nwin, lwin, linc);
    
    cout << "Semblance is " << semb << '\n';

    return 0;
}
