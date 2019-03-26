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

#include <boost/timer.hpp>
#include <iostream>
#include <list>

#include <boost/numeric/linear_algebra/accumulate.hpp>
#include <boost/numeric/linear_algebra/concept_maps.hpp>
#include <boost/numeric/linear_algebra/operators.hpp>


template <typename Element>
void test_accumulate(const char* name)
{
    const int   array_size= 10;
    Element     array[array_size];
    for (int i= 0; i < array_size; i++) 
    	array[i]= Element(i);

    std::list<Element> l;
    for (int i= 0; i < array_size; i++) 
    	l.push_back(Element(i));
    
    std::cout << '\n' << name << '\n' << " Add: ";
    math::accumulate(&array[0], array+array_size, Element(0), math::add<Element>());
    std::cout << "Mult: ";
    math::accumulate(array, array+array_size, Element(1), math::mult<Element>());
    std::cout << "Mult [with a list]: ";
    math::accumulate(l.begin(), l.end(), Element(1), math::mult<Element>());
    std::cout << " Min: ";
    math::accumulate(array, array+array_size, Element(1000), math::min<Element>());
    std::cout << " Max: ";
    math::accumulate(array, array+array_size, Element(-1000), math::max<Element>());
}


int main(int, char* [])
{
    test_accumulate<int>("int");
    test_accumulate<float>("float");
    test_accumulate<double>("double");
    std::cout << '\n';

    return 0;
}
