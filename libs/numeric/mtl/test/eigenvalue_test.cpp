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
//
// Written by Marc Hartung

#include <boost/numeric/mtl/mtl.hpp>
#include <cmath>


int main() {
    mtl::mat::dense2D<double> A(10, 10);
    
    int ev_count = 0;
    
    // The eigenvalues are {10, 20, 30, ..., 100}
    
    A = 54.172, -14.5819, -7.74761, -7.95597, 1.60174, -17.9864, 3.27844, 11.7596, -1.28794, 1.87001, 
	-11.6764, 53.1245, -4.85491, 4.32814, -4.33267, 1.43186, 28.9238, 6.22295, -10.7843, 3.84623, 
	-10.0254, -4.20435, 44.4108, -19.8427, -9.92564, 1.47404, -4.08042, -11.4571, 3.06065, -7.1122, 
	-12.14, 1.97065, -17.368, 75.8573, 3.3052, -0.552501, 0.376509, -1.99455, 14.517, -8.8973, 
	-1.62797, -6.57477, -8.85886, 3.76852, 48.0681, 3.16219, -4.7653, 9.94849, 11.9553, 1.03926, 
	-18.4308, -2.98651, 0.548431, 0.0521285, 9.22361, 61.9968, -7.21312, 9.29003, 0.117305, 1.61741, 
	3.63301, 26.2867, -4.06843, 2.81207, -7.74316, -11.215, 51.3995, -4.16883, -7.48933, 2.27791, 
	7.82304, 3.30253, -13.8787, 2.12292, 10.5754, 6.73106, -4.2492, 41.0607, 18.4729, -6.49337, 
	-3.39889, -10.3545, -0.828751, 17.1083, 6.65637, -1.54847, -9.39066, 19.2665, 59.9372, 8.98516, 
	0.991482, 2.81724, -1.99112, -10.9052, 1.77954, 3.20488, 0.519006, -4.37109, 9.20029, 59.9732;
	
    mtl::mat::eigenvalue_solver<mtl::mat::dense2D<double> > ES(A);
    ES.calc();
    for(int i=1;i<=10;i++) {
	for(int j=0;j<10;j++) {
	    if(std::abs(ES.get_eigenvalues()[j].real()-(i*10)) < 1.0e-4) {
		ev_count++;
		break;
	    }
	}
    }
    
    MTL_THROW_IF(ev_count != 10, mtl::runtime_error("Should be 10 eigenvalues."));
        
    return 0;
}
