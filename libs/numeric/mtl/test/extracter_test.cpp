#include <boost/numeric/mtl/mtl.hpp>
#include <stdexcept>
#include <sstream>

using namespace mtl;

template< typename T >
void check(const T& des, T g) {
	if(des != g) {
		std::stringstream ss;
		ss<<" got: "<<des<<"\n wanted: "<<g<<"\n ERROR";
 		//std::cout<< "Error= " << ss.str() << "\n";
		throw std::runtime_error(ss.str());
	}
}

template< typename Vector >
void test(const Vector& vec) {
	typedef typename Vector::value_type value_type;
	value_type one, two, four;
	{
		mtl::vec::extracter< Vector > extract(vec);
		extract[1] >> one;
		extract[2] >> two;
		extract[3] >> four;
	}
	check(one, 1.0);
	check(two, 2.0);
	check(four, 4.0);
}

template< typename Vector >
void test(const Vector& , typename mtl::vec::extracter< Vector >::buffer_type& buffer) {
	typedef typename Vector::value_type value_type;
	value_type one, two, four;
	{
		mtl::vec::extracter< Vector > extract(buffer);
		extract[1] >> one;
		extract[2] >> two;
		extract[3] >> four;
	}
	check(one, 1.0);
	check(two, 2.0);
	check(four, 4.0);
}

int main(int , char** ) {

	dense_vector< double > dense(4);
	dense[1]=1.0;
	dense[2]=2.0;
	dense[3]=4.0;

	test(dense);
	mtl::vec::extracter< dense_vector< double > >::buffer_type buffer(dense);
	test(dense, buffer);

	return 0;
}

