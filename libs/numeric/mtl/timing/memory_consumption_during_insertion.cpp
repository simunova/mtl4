#include <iostream>
#include <vector>
#include <unistd.h>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/timer.hpp>

using namespace std;

long bw(long w) 
{ 
    long k= 0;
    //for (long long i0= 0; i0 < w; i0++)
      for (long long i= 0; i < w; i++)
	for (long j= 0; j < 100000000; j++)
	    k+= j;
      std::cout << k << endl;
      return k;
}

int main()
{
    const int s= 1000000, n= 4, z= 32, t= 15;
    long w= 100, r;


#if 0
    std::vector<double> v(10000000);
    std::cout << "Vector 10M created. Capacity == " << v.capacity() << endl;
    r= bw(w);

    v.resize(20000000);
    std::cout << "Vector auf 20M vergroessert. Capacity == " << v.capacity() << endl;
    r= bw(w);

    v.resize(5000000);
    std::cout << "Vector auf 5M verkleinert. Capacity == " << v.capacity() << endl;
    r= bw(w);

    {
	std::vector<double> w(v.begin(), v.end());
	swap(v, w);
    }
    std::cout << "Mit 5M geswappt. Capacity == " << v.capacity() << endl;
    r= bw(w);
#endif

    

    mtl::compressed2D<double> A(s, s);
    std::cout << "Matrix created." << endl;
    r= bw(w);

    {
	mtl::mat::inserter<mtl::compressed2D<double> > ins(A, z);
	std::cout << "Inserter created." << endl;
	r= bw(w);
 
	std::cout << "start insertion." << endl;
	for (int i= n; i < s; i++)
	    for (int j= i - n; j < i; j++)
		ins[i][j] << 9.0;
	std::cout << "Matrix filled." << endl;
	r= bw(w);
    }
    std::cout << "Instructor destroyed." << endl;
    r= bw(w);

    A.shrink();
    std::cout << "After A.shrink()." << endl;
    r= bw(w);

    return 0;
}
