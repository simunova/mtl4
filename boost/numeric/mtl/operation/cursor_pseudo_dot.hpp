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

#ifndef MTL_CURSOR_PSEUDO_DOT_INCLUDE
#define MTL_CURSOR_PSEUDO_DOT_INCLUDE

namespace mtl {



namespace functor {

    template <unsigned MaxDepth, typename Value, typename Cursor1, typename Prop1, typename Cursor2, typename Prop2, unsigned Depth>
    struct cursor_pseudo_dot_block
    {
	static unsigned const offset= MaxDepth - Depth;
	
	void operator() (Cursor1 i1, Prop1& prop1, Cursor2 i2, Prop2& prop2,
			 Value& s0, Value& s1, Value& s2, Value& s3,
			 Value& s4, Value& s5, Value& s6, Value& s7)
	{
	    Cursor1 tmp1(i1); tmp1+= offset;
	    Cursor2 tmp2(i2); tmp2+= offset;
	    s0+= prop1(*tmp1) * prop2(*tmp2);
	    // s0+= prop1(i1 + offset) * prop2(i2 + offset);
	    typedef cursor_pseudo_dot_block<MaxDepth, Value, Cursor1, Prop1, Cursor2, Prop2, Depth-1> block_rest;
	    block_rest() (i1, prop1, i2, prop2, s1, s2, s3, s4, s5, s6, s7, s0);
	}
    };

    //template <>
    template <unsigned MaxDepth, typename Value, typename Cursor1, typename Prop1, typename Cursor2, typename Prop2>
    struct cursor_pseudo_dot_block<MaxDepth, Value, Cursor1, Prop1, Cursor2, Prop2, 1>
    {
	static unsigned const offset= MaxDepth - 1;
	
	void operator() (Cursor1 i1, Prop1& prop1, Cursor2 i2, Prop2& prop2,
			 Value& s0, Value&, Value&, Value&,
			 Value&, Value&, Value&, Value&)
	{
	    s0+= prop1(*(i1 + offset)) * prop2(*(i2 + offset));
	}
    };      

    template <unsigned MaxDepth, typename Value, typename Cursor1, typename Prop1, typename Cursor2, typename Prop2>
    struct cursor_pseudo_dot_t
    {
	Value operator() (Cursor1 i1, Cursor1 end1, Prop1& prop1, Cursor2 i2, Prop2& prop2)
	{
	    using math::zero;
	    Value         ref, my_zero(zero(ref)),
                          s0= my_zero, s1= my_zero, s2= my_zero, s3= my_zero, 
		          s4= my_zero, s5= my_zero, s6= my_zero, s7= my_zero;
	    std::size_t size= end1 - i1, blocks= size / MaxDepth, blocked_size= blocks * MaxDepth;

	    typedef cursor_pseudo_dot_block<MaxDepth, Value, Cursor1, Prop1, Cursor2, Prop2, MaxDepth> dot_block_type;
	    for (unsigned i= 0; i < blocked_size; i+= MaxDepth, i1+= MaxDepth, i2+= MaxDepth) {
		dot_block_type()(i1, prop1, i2, prop2, s0, s1, s2, s3, s4, s5, s6, s7);
	    }

	    typedef cursor_pseudo_dot_block<MaxDepth, Value, Cursor1, Prop1, Cursor2, Prop2, MaxDepth> dot_single_type;
	    s0+= s1 + s2 + s3 + s4 + s5 + s6 + s7;
	    for (unsigned i= blocked_size; i < size; ++i, ++i1, ++i2)
		dot_single_type()(i1, prop1, i2, prop2, s0, s1, s2, s3, s4, s5, s6, s7);
	    return s0;
	}
    };
 
} // namespace functor 

template <unsigned MaxDepth, typename Value, typename Cursor1, typename Prop1, typename Cursor2, typename Prop2>
Value cursor_pseudo_dot(Cursor1 i1, Cursor1 end1, Prop1 prop1, Cursor2 i2, Prop2 prop2, Value)
{
    return functor::cursor_pseudo_dot_t<MaxDepth, Value, Cursor1, Prop1, Cursor2, Prop2>()(i1, end1, prop1, i2, prop2);
}



} // namespace mtl

#endif // MTL_CURSOR_PSEUDO_DOT_INCLUDE
