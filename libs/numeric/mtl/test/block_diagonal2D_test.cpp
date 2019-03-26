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

using namespace mtl;

static const int DIM =3;
static const int NB  =4;
typedef dense2D<double, mat::parameters<tag::row_major,
					   mtl::index::c_index, mtl::fixed::dimensions< DIM, DIM> > > Mat3;
typedef mat::block_diagonal2D<Mat3> MatB;

int main()
{
    double array[][3]={{1.0, 0.0, 5.5}, {0.0, 2.0, 4.0}, {0.0, 0.0, 3.0}};
    Mat3 E3(array);

    MatB Eb(NB*DIM, NB*DIM);

    for(int i=0; i<NB; i++)
	Eb.insert(i*DIM, (i+1)*DIM, E3) ;

    mtl::io::tout << "E3 = \n" << E3 << "\n";  
    mtl::io::tout << "Eb = \n" << Eb << "\n";   

    MTL_THROW_IF((Eb(9, 9) != 1.0), mtl::unexpected_result());
    MTL_THROW_IF((Eb(9, 10) != 0.0), mtl::unexpected_result());
    
    return 0;
}

