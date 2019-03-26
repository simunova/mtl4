// Filename: insert_scope.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

class world_matrix
{
    typedef mat::inserter<compressed2D<double>, update_plus<double> > inserter_type;
  public:
    world_matrix(unsigned nrows, unsigned ncols) : A(nrows, ncols), ins(0) {}

    void add_entry(unsigned row, unsigned col, double value) 
    {  (*ins)[row][col] << value;    }

    void start_insertion() // must be called before first insertion
    { if (!ins) ins= new inserter_type(A); }

    void finish_insertion() // must be called before first usage
    {
	if (ins) {
	    delete ins;
	    ins= 0;
	}
    }

    friend inline std::ostream& operator<<(std::ostream& os, const world_matrix& w)
    { return os << w.A; }

  private:
    compressed2D<double>      A;
    inserter_type*            ins;
};


int main(int, char**)
{
    world_matrix              A(3, 3);

    A.start_insertion();
    A.add_entry(0, 0, 2.0);
    A.add_entry(1, 2, 0.5);
    A.add_entry(2, 1, 3.0);
    A.finish_insertion();

    std::cout << "A is \n" << A;  // we can access A now

    return 0;
}

