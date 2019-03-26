// Filename: vampir_example.cpp (part of MTL4)

#include <boost/numeric/mtl/mtl.hpp>

using namespace std;  

template <typename Vector>
inline void my_add(const Vector& u, const Vector& v, Vector& w)
{
    mtl::vampir_trace<12800> tracer;
    w= u + v;
}

namespace mtl { namespace vpt { 
    template <> std::string vampir_trace<12800>::name("my_add");
}}

int main(int, char**)
{
    mtl::vampir_trace<9999> tracer;
    mtl::dense_vector<float> v1(3, 1.0), v2(3, 2.0), v3;
    my_add(v1, v2, v3);

    return 0;
}
