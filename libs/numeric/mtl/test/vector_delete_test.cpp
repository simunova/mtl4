// MTL4 Test 05.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

typedef mtl::dense_vector<double> mtlVec;
mtlVec testFunc( mtlVec A,double Scale);

int main(int, char**)
{
    mtlVec u(3);
    u(0)=1.1;
    u(1)=2.2;
    u(2)=3.3;
    mtlVec myVec(3);
    myVec=testFunc(u,55.5);
    print_vector(myVec);

    return 0;
}


mtlVec testFunc(mtlVec A ,double Scale)
{
    const std::size_t testSize=size(A);
    mtlVec testVec(testSize);
    testVec=Scale*A;
    print_vector(testVec);

    return testVec;
}
