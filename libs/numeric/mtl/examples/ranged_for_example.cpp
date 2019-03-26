// Filename: ranged_for_example.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main()
{
#if defined(MTL_WITH_RANGEDFOR) && defined(MTL_WITH_INITLIST)
    using mtl::irange;
    for (int i : irange(3, 9))
	std::cout << i << '\n';

    mtl::dense_vector<double> v= {3, 4, 5, 7};
    for (int i : irange(size(v)))
	v[i]+= 2.5;
    std::cout << "v is " << v << '\n';
#endif

    return 0;
}
