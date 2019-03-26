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

#ifndef MTL_MAKE_COPY_OR_REFERENCE_INCLUDE
#define MTL_MAKE_COPY_OR_REFERENCE_INCLUDE

namespace mtl {


/// Helper class to avoid avoidable copies for input parameters
/** Container is referred if it has already target type, otherwise copied.
    Create an object of this type and pass the value member variable to the function,
    e.g. make_in_copy_or_reference<Tgt, Src> copy_or_ref(v); f(copy_or_ref.value);
    where Src is the type of v and Tgt the type of f's argument.
**/
template <typename Target, typename Source>
struct make_in_copy_or_reference
{
    typedef Target             type;
    explicit make_in_copy_or_reference(const Source& src) : value(src) {}
    Target value;
};

template <typename Target>
struct make_in_copy_or_reference<Target, Target>
{
    typedef const Target&      type;
    explicit make_in_copy_or_reference(const Target& src) : value(src) {}
    const Target& value;
};


/// Helper class to avoid avoidable copies for output parameters.
/** Container is referred if it has already target type, otherwise copied at destruction.
    Create an object of this type and pass the value member variable to the function,
    e.g. make_in_copy_or_reference<Tgt, Src> copy_or_ref(v); f(copy_or_ref.value);
    where Src is the type of v and Tgt the type of f's argument.
    Target must be DefaultConstructible.
**/
template <typename Target, typename Source>
struct make_out_copy_or_reference
{
    explicit make_out_copy_or_reference(Source& src) : src(src) {}
    ~make_out_copy_or_reference() { src= value; }

    Target  value;
private:
    Source& src;
};

template <typename Target>
struct make_out_copy_or_reference<Target, Target>
{
    explicit make_out_copy_or_reference(Target& src) : value(src) {}
    Target& value;
};


/// Helper class to avoid avoidable copies for input-output parameters.
/** Container is referred if it has already target type, otherwise copied construction and destruction.
    Create an object of this type and pass the value member variable to the function,
    e.g. make_in_copy_or_reference<Tgt, Src> copy_or_ref(v); f(copy_or_ref.value);
    where Src is the type of v and Tgt the type of f's argument.
**/
template <typename Target, typename Source>
struct make_in_out_copy_or_reference
{
    explicit make_in_out_copy_or_reference(Source& src) : value(src), src(src) {}
    ~make_in_out_copy_or_reference() { src= value; }

    Target  value;
private:
    Source& src;
};

template <typename Target>
struct make_in_out_copy_or_reference<Target, Target>
{
    explicit make_in_out_copy_or_reference(Target& src) : value(src) {}
    Target& value;
};

} // namespace mtl

#endif // MTL_MAKE_COPY_OR_REFERENCE_INCLUDE
