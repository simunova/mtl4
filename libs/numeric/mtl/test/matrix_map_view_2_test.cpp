/*
 *  matrix_map_view_test_2.cpp
 *  MTL
 *
 *	Test mat::rscaled_view and mat::divide_by_view
 *
 *  Created by Hui Li (huil@Princeton.EDU)
 *
 */

#include <iostream>
#include <cmath>
#include <complex>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp> 
#include <boost/numeric/mtl/matrix/compressed2D.hpp> 
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/matrix/hermitian_view.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/recursion/predefined_masks.hpp>
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/mtl/operation/scale.hpp>
#include <boost/numeric/mtl/operation/rscale.hpp>
#include <boost/numeric/mtl/operation/divide_by.hpp>
#include <boost/numeric/mtl/operation/hermitian.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>
#include <boost/numeric/mtl/operation/mult_result.hpp>
#include <boost/numeric/mtl/operation/mult.hpp>
#include <boost/numeric/mtl/operation/exp.hpp>
#include <boost/numeric/mtl/operation/iota.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>


using namespace std;  

typedef complex<double> ct;

double value(double)
{
    return 7.0;
}

complex<double> value(complex<double>)
{
    return ct(7.0, 1.0);
}

// right scaled value
double rsvalue(double)
{
    return 14.0;
}

ct rsvalue(ct)
{
    return ct(14.0, 2.0);
}

// complex scaled value
ct crsvalue(double)
{
    return ct(0.0, 7.0);
}

ct crsvalue(ct)
{
    return ct(-1.0, 7.0);
}


template <typename Matrix>
void test(Matrix& matrix, const char* name)
{
    set_to_zero(matrix);
    typedef typename Matrix::value_type value_type;
    value_type ref(0);
	
    {
		mtl::mat::inserter<Matrix>  ins(matrix);
		ins(2, 3) << value(ref);
		ins(4, 3) << value(ref) + 1.0;
		ins(2, 5) << value(ref) + 2.0;
    }
	
    cout << "\n\n" << name << "\n";
    cout << "Original matrix:\n" << matrix << "\n";
	
#if 0
    std::cout << "ashape double: " << typeid(typename ashape::ashape<double>::type).name() << "\n";
    std::cout << "ashape matrix: " << typeid(typename ashape::ashape<Matrix>::type).name() << "\n";
    typedef typename ashape::ashape<double>::type   dshape;
    typedef typename ashape::ashape<Matrix>::type   mshape;
	
    std::cout << "ashape mult_op: " << typeid(typename ashape::mult_op<dshape, mshape>::type()).name() << "\n";
    typedef typename ashape::mult_op<dshape, mshape>::type rshape;
	
    
    std::cout << "result type: " << typeid(typename traits::mult_result_aux<double, Matrix, rshape>::type).name() << "\n";
#endif
	
	
	// test rscaled_view
    mtl::mat::rscaled_view<Matrix, double>  rscaled_matrix(matrix, 2.0);
    cout << "matrix  right scaled with 2.0\n" << rscaled_matrix << "\n";
    MTL_THROW_IF(rscaled_matrix(2, 3) != rsvalue(ref), mtl::runtime_error("right scaling wrong"));
	
    cout << "matrix  right scaled with 2.0 (as operator)\n" << matrix * 2.0 << "\n";
    MTL_THROW_IF((matrix * 2.0)(2, 3) != rsvalue(ref), mtl::runtime_error("right scaling wrong"));
	
	
    mtl::mat::rscaled_view<Matrix, ct>  crscaled_matrix(matrix, ct(0.0, 1.0));
    cout << "matrix right scaled with i (complex(0, 1))\n" << crscaled_matrix << "\n";
    MTL_THROW_IF(crscaled_matrix(2, 3) != crsvalue(ref), mtl::runtime_error("complex right scaling wrong"));
	
    cout << "matrix right scaled with 2.0 (free function)\n" << rscale(matrix, 2.0) << "\n";
    MTL_THROW_IF(rscale(matrix, 2.0)(2, 3) != rsvalue(ref), mtl::runtime_error("scaling wrong"));
	
    cout << "matrix right scaled with i (complex(0, 1)) (free function)\n" << rscale(matrix, ct(0.0, 1.0)) << "\n";
    MTL_THROW_IF(rscale(matrix, ct(0.0, 1.0))(2, 3) != crsvalue(ref), mtl::runtime_error("complex right scaling wrong"));
	
	// test divide_by_view
    mtl::mat::divide_by_view<Matrix, double>  div_matrix(matrix, 0.5);
    cout << "matrix divide by 0.5\n" << div_matrix << "\n";
    MTL_THROW_IF(div_matrix(2, 3) != rsvalue(ref), mtl::runtime_error("divide_by wrong"));
	
    cout << "matrix divide by 0.5 (as operator)\n" << matrix / 0.5 << "\n";
    MTL_THROW_IF((matrix / 0.5)(2, 3) != rsvalue(ref), mtl::runtime_error("divide_by wrong"));
	
	
    mtl::mat::divide_by_view<Matrix, ct>  cdiv_matrix(matrix, ct(0.0, -1.0));
    cout << "matrix divide by -i (complex(0, -1))\n" << cdiv_matrix << "\n";
    MTL_THROW_IF(cdiv_matrix(2, 3) != crsvalue(ref), mtl::runtime_error("complex divide_by wrong"));
	
    cout << "matrix divide by 0.5 (free function)\n" << rscale(matrix, 0.5) << "\n";
    MTL_THROW_IF(divide_by(matrix, 0.5)(2, 3) != rsvalue(ref), mtl::runtime_error("scaling wrong"));
	
    cout << "matrix divide by -i (complex(0, 1)) (free function)\n" << divide_by(matrix, ct(0.0, -1.0)) << "\n";
    MTL_THROW_IF(divide_by(matrix, ct(0.0, -1.0))(2, 3) != crsvalue(ref), mtl::runtime_error("complex right scaling wrong"));

    mtl::mat::exp_view<Matrix>   exp_matrix(matrix);
    // cout << "exp matrix is\n" << exp_matrix << '\n';
    cout << "exp matrix is\n" << exp(matrix) << '\n';
    
    mtl::dense_vector<value_type> v(num_cols(matrix), value_type(1)), w;
    
    //  To be deleted some day
//     cout <<  typeid(typename mtl::sfunctor::exp<typename Matrix::value_type>::result_type).name() << '\n';
//     cout <<  typeid(typename mtl::mat::map_view<mtl::sfunctor::exp<typename Matrix::value_type>, Matrix>::value_type).name() << '\n';
//     cout <<  typeid(typename mtl::mat::exp_view<Matrix>::value_type).name() << '\n';
//     cout <<  typeid(typename mtl::Collection<mtl::mat::exp_view<Matrix> >::value_type).name() << '\n';
//     cout <<  typeid(typename mtl::Collection<mtl::mat::imag_view<Matrix> >::value_type).name() << '\n';
    
    w= exp(matrix) * v;
    cout << "exp(A) * vec(1) = " << w << '\n';
    // cout << "exp(A) * vec(1) = " << exp(matrix) * v << '\n'; // maybe later

}



int main(int argc, char* argv[])
{
    using namespace mtl;
    unsigned size= 7; 
    if (argc > 1) size= atoi(argv[1]); 
	
    dense2D<double>                                      dr(size, size);
    dense2D<double, mat::parameters<col_major> >      dc(size, size);
    morton_dense<double, recursion::morton_z_mask>       mzd(size, size);
    morton_dense<double, recursion::doppled_2_row_mask>  d2r(size, size);
    compressed2D<double>                                 cr(size, size);
    compressed2D<double, mat::parameters<col_major> > cc(size, size);
	
    dense2D<complex<double> >                            drc(size, size);
    compressed2D<complex<double> >                       crc(size, size);
	
	
    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
    test(cr, "Compressed row major");
    test(cc, "Compressed column major");
    test(drc, "Dense row major complex");
    test(crc, "Compressed row major complex");
	
    return 0;
}
