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

#ifndef MTL_IO_MATRIX_FILE_INCLUDE
#define MTL_IO_MATRIX_FILE_INCLUDE

namespace mtl { namespace io {

template <typename MatrixIFStream, typename MatrixOFStream>
class matrix_file
{
  public:
    explicit matrix_file(const std::string& fname) : fname(fname) {}
    explicit matrix_file(const char* fname) : fname(fname) {}

    std::string file_name() const { return fname; }

    template <typename Collection>
    matrix_file& operator=(const Collection& c)
    {
	MatrixOFStream stream(fname);
	stream << c;
	return *this;
    }

  protected:
    std::string fname;
};

}} // namespace mtl::io

#endif // MTL_IO_MATRIX_FILE_INCLUDE
