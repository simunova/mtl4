#include <iostream>
#include <cmath>
#include <vector>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;
    dense_vector<double>                                        cd(5, 1.0);
    std::vector<dense_vector<double> >                          v(5, cd);
 
    for (unsigned i= 0, c= 1; i < size(v); ++i)
	for (unsigned j= 0; j < size(v[i]); ++j, c++)
	    v[i][j]= double((i + j) % 5);

    std::cout << "w initially\n";
    std::vector<dense_vector<double> > 				w(v),x(v);
    for (unsigned i= 0; i < size(w); ++i)
	std::cout << w[i] << "\n";
 
    orth(w);
    std::cout << "\nw orthogonalized\n";
    for (unsigned i= 0; i < size(w); ++i)
	std::cout << w[i] << "\n";
    
    std::cout << "\nTest: dot product of orthogonalized w\n";
    for (unsigned i= 0, c= 1; i < size(w); ++i) {
	for (unsigned j= 0; j < size(w); ++j, ++c)
	    std::cout << (dot(w[i], w[j])<1.e-8 ? 0 : dot(w[i], w[j]) )<< " " ;
	std::cout << "\n";
    }   
    std::cout << "\nThe according factors are: \n" << orthogonalize_factors(v) << '\n';
    
    orth(x,0); orth(x,1);
    std::cout << "\nw(0), w(1) orthogonalized\n";
    for (unsigned i= 0; i < size(w); ++i)
	std::cout << x[i] << "\n";
    

    return 0;
}
