// Filename: insert_class_expensive.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

class world_matrix
{
  public:
    world_matrix(unsigned nrows, unsigned ncols) : A(nrows, ncols) {}

    void add_entry(unsigned row, unsigned col, double value) 
    {
	// Extremely expensive -> must not be done 
	mat::inserter<compressed2D<double>, update_plus<double> > ins(A); 
	ins[row][col] << value;
    }

    friend inline std::ostream& operator<<(std::ostream& os, const world_matrix& w)
    { return os << w.A; }

  private:
    compressed2D<double>              A;
};


int main(int, char**)
{
    world_matrix              A(3, 3);

    A.add_entry(0, 0, 2.0);
    A.add_entry(1, 2, 0.5);
    A.add_entry(2, 1, 3.0);

    std::cout << "A is \n" << A;  // we can access A now

    return 0;
}

