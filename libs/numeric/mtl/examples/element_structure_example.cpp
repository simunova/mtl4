#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;
int main(int, char**)
{
    typedef double value_type;
    typedef int    size_type;
    const int nb_elements= 6, nb_nodes= 12;
    
    value_type array[][4]= {{2, 3,   4,   5}, 
			    {4, 10, 13,  16},
			    {6, 25, 38,  46},
			    {8, 32, 77, 100}};
    mtl::dense2D<value_type>   E_mat(array);
    
    std::cout<<"E_mat=\n"<< E_mat <<"\n";

    mtl::mat::element_structure<value_type> A;
    
    typedef mtl::mat::element<value_type>	element_type;
    element_type* elements = new element_type[nb_elements];
    
    mtl::dense_vector<size_type>  index_a(4, 0), 
				   index_b(4, 0), 
				   index_c(4, 0), 
				   index_d(4, 0), 
				   index_e(4, 0), 
				   index_f(4, 0);
	
    // construct nodes for every element
    index_a[0]= 0; index_a[1]= 1; index_a[2]= 4; index_a[3]= 5;
    index_b= index_a + 1;
    index_c= index_a + 2;
    index_d= index_a + 4;
    index_e= index_a + 5;
    index_f= index_a + 6;
   
    //construct the 6 elements from the example grid
    element_type a(0, index_a, E_mat);
    element_type b(1, index_b, E_mat);
    element_type c(2, index_c, E_mat);
    element_type d(3, index_d, E_mat);
    element_type e(4, index_e, E_mat);
    element_type f(5, index_f, E_mat);
    
    //construct neighborhood information for each element
    a.add_neighbors(&b, &d, &e);
    b.add_neighbors(&a, &c, &d, &e, &f);
    c.add_neighbors(&a, &b, &e);
    d.add_neighbors(&a, &b, &e);
    e.add_neighbors(&a, &b, &c, &d, &f);
    f.add_neighbors(&b, &c, &e);

    std::cout<< "a=" << a << "\n";
    
    //construct array of elements
    elements[0]=a;
    elements[1]=b;
    elements[2]=c;
    elements[3]=d;
    elements[4]=e;
    elements[5]=f;
    
    //construct element_structure from the 6 single elements
    A.consume(nb_elements, nb_nodes, elements);
    mtl::dense_vector<value_type> x(nb_nodes, 1.0), test(A * x);
    
    std::cout<< "test="<< test << "\n";
    
    return 0;
}
