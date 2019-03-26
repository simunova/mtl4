// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>


using namespace std;  


namespace std {

    struct stream
    {
	stream(std::ostream &s= std::cout) : s(s) {}
	
	template <typename T> 
	stream& operator<<(const T& x) { s << "Ciao " << x; }
	
	void flush() { s.flush(); }

	std::ostream &s;
    };
    
    inline stream& endl(stream& s) { s << "\n"; s.flush(); }
}

void f() {}

int main(int argc, char* argv[])
{
    std::stream my_stream;
    my_stream << "bella " << 79;
    my_stream << f;
    //my_stream << std::endl;

    return 0;
}
