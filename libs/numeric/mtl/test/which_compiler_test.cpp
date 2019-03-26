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


int main(int, char**)
{
# if defined(__INTEL_COMPILER)
    std::cout << "Intel compiler, version (intern) " << __INTEL_COMPILER << ", i.e. icc/icpc " 
	      << double(__INTEL_COMPILER) / 100.0 << '\n';
# elif defined(__clang__)
    std::cout << "Clang C++ compiler, version " << __clang_major__ << '.' << __clang_minor__ << '\n';
# elif defined(__GNUC__)
    std::cout << "GNU compiler, i.e. g++ " << __GNUC__ << '.' << __GNUC_MINOR__ << '.' << __GNUC_PATCHLEVEL__ << '\n';
# elif defined(__PGI)
    std::cout << "Portland group compiler (MTL4 not tested for it yet).\n";
# elif defined(__IBMCPP__)
    std::cout << "IBM compiler xlc++ (MTL4 not tested for it yet).\n"; 
    // version __xlC__ (or __xlc__) VVRM version, release, modification
# elif defined(_MSC_VER)
    std::cout << "Visual studio C++ compiler, version (intern) " << _MSC_VER << '\n';
# else
    std::cout << "Unknown compiler\n";
# endif

#if 0 // pgCC does not terminate with this test, Intel compiler ambiguous
    int compilers= 0;
# ifdef __INTEL_COMPILER
    compilers++;
# endif
# ifdef __GNUC__
    compilers++;
# endif
# ifdef __PGI
    compilers++;
# endif
# ifdef __IBMCPP__
    compilers++;
# endif
# ifdef _MSC_VER
    compilers++;
# endif
# ifdef __clang__
    compilers++;
# endif
    if (compilers > 1) 
	std::cout << "More then one compiler detected!!!!\n";
#endif 

    return 0;
}
