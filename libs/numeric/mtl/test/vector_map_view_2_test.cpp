/*
 *  vector_map_view_test_2.cpp
 *  MTL
 *
 *	Test rscaled_view and divide_by_view
 *
 *  Created by Hui Li (huil@Princeton.EDU)
 *
 */


#include <iostream>
#include <cmath>
#include <complex>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/mtl/operation/rscale.hpp>
#include <boost/numeric/mtl/operation/divide_by.hpp>

#if 0
#include <boost/numeric/mtl/operation/hermitian.hpp>
#endif


using std::cout;  using std::complex;

typedef complex<double> ct;

double value(double)
{
    return 7.0;
}

complex<double> value(complex<double>)
{
    return ct(7.0, 1.0);
}

// rscaled value
double rsvalue(double)
{
    return 14.0;
}

ct rsvalue(ct)
{
    return ct(14.0, 2.0);
}

// complex rscaled value
ct crsvalue(double)
{
    return ct(0.0, 7.0);
}

ct crsvalue(ct)
{
    return ct(-1.0, 7.0);
}


template <typename Vector>
void test(Vector& vector, const char* name)
{
    set_to_zero(vector);
    typename Vector::value_type ref(0);
	
#if 1	
    vector[2]= value(ref);
    vector[4]= value(ref) + 1.0;
    vector[5]= value(ref) + 2.0;
	
#else // When sparse vectors are used there should be an inserter class for vectors too
    {
		inserter<Vector>  ins(vector);
		ins(2) << value(ref);
		ins(4) << value(ref) + 1.0;
		ins(5) << value(ref) + 2.0;
    }
#endif
	
    cout << "\n\n" << name << "\n";
    cout << "Original vector:\n" << vector << "\n";
	
	// test rscaled_view
    mtl::vec::rscaled_view<Vector,double>  rscaled_vector(vector,2.0);
    cout << "vector right scaled with 2.0\n" << rscaled_vector << "\n";
    MTL_THROW_IF(rscaled_vector(2) != rsvalue(ref), mtl::runtime_error("right scaling wrong"));
    
    mtl::vec::rscaled_view<Vector,ct>  crscaled_vector(vector,ct(0.0, 1.0));
    cout << "vector right scaled with i (complex(0, 1))\n" << crscaled_vector << "\n";
    MTL_THROW_IF(crscaled_vector(2) != crsvalue(ref), mtl::runtime_error("complex right scaling wrong"));
	
    cout << "vector right scaled with 2.0 (free function)\n" << rscale(vector,2.0) << "\n";
    MTL_THROW_IF(rscale(vector,2.0)(2) != rsvalue(ref), mtl::runtime_error("right scaling wrong"));
	
    cout << "vector right scaled with i (complex(0, 1)) (free function)\n" << rscale(vector,ct(0.0, 1.0)) << "\n";
    MTL_THROW_IF(rscale(vector,ct(0.0, 1.0))(2) != crsvalue(ref), mtl::runtime_error("complex right scaling wrong"));
	
	// test divide_by_view
    mtl::vec::divide_by_view<Vector,double>  div_vector(vector,0.5);
    cout << "vector divide by 0.5\n" << div_vector << "\n";
    MTL_THROW_IF(div_vector(2) != rsvalue(ref), mtl::runtime_error("divide_by wrong"));
    
    mtl::vec::divide_by_view<Vector,ct>  cdiv_vector(vector,ct(0.0, -1.0));
    cout << "vector divide by -i (complex(0, -1))\n" << cdiv_vector << "\n";
    MTL_THROW_IF(cdiv_vector(2) != crsvalue(ref), mtl::runtime_error("complex divide_by wrong"));
	
    cout << "vector divide by 0.5 (free function)\n" << divide_by(vector,0.5) << "\n";
    MTL_THROW_IF(divide_by(vector,0.5)(2) != rsvalue(ref), mtl::runtime_error("divide_by wrong"));
	
    cout << "vector divide by -i (complex(0, -1)) (free function)\n" << divide_by(vector,ct(0.0, -1.0)) << "\n";
    MTL_THROW_IF(divide_by(vector,ct(0.0, -1.0))(2) != crsvalue(ref), mtl::runtime_error("complex divide_by wrong"));
	
}



int main(int argc, char* argv[])
{
    unsigned size= 7; 
    if (argc > 1) size= atoi(argv[1]); 
	
    mtl::dense_vector<double>                                 dv(size);
    mtl::dense_vector<complex<double> >                       drc(size);
	
    test(dv, "Dense double vector");
    test(drc, "Dense complex vector");
	
    return 0;
}
