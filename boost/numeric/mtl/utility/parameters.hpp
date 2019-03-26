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

#ifndef MTL_TRAITS_PARAMETERS_INCLUDE
#define MTL_TRAITS_PARAMETERS_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>

namespace mtl { namespace traits {

template <typename T>
struct parameters
{
    typedef typename T::parameters  type;
};

template <typename Vector>
struct parameters<mtl::mat::multi_vector<Vector> >
{
    typedef mtl::mat::parameters<>      type;
};

template <typename Functor>
struct parameters<mtl::mat::implicit_dense<Functor> >
{
    typedef mtl::mat::parameters<>      type;
};

template <typename Value>
struct parameters<mtl::mat::ones_matrix<Value> >
  : public parameters<mtl::mat::implicit_dense<mtl::mat::ones_functor<Value> > > 
{};

template <typename Value>
struct parameters<mtl::mat::hilbert_matrix<Value> >
  : public parameters<mtl::mat::implicit_dense<mtl::mat::hilbert_functor<Value> > > 
{};

template <typename Vector1, typename Vector2>
struct parameters<mtl::mat::outer_product_matrix<Vector1, Vector2> >
  : public parameters<mtl::mat::implicit_dense<mtl::mat::outer_product_functor<Vector1, Vector2> > > 
{};


}} // namespace mtl::traits

#endif // MTL_TRAITS_PARAMETERS_INCLUDE
