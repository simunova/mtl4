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

#ifndef MTL_MULTI_ACTION_BLOCK_INCLUDE
#define MTL_MULTI_ACTION_BLOCK_INCLUDE

namespace mtl {

/* 
     Functor for unrolling arbitrary loops with independent operations
     ...

 */
template <typename MultiAction, unsigned Steps> struct multi_action_block;


template <typename MultiAction, unsigned MaxSteps, unsigned RemainingSteps>
struct multi_action_helper
{
    static unsigned const step= MaxSteps - RemainingSteps;

    void operator() (MultiAction const& action) const
    {
	action(step);
	multi_action_helper<MultiAction, MaxSteps, RemainingSteps-1>()(action);
    }

    void operator() (MultiAction& action) const
    {
	action(step);
	multi_action_helper<MultiAction, MaxSteps, RemainingSteps-1>()(action);
    }    
};


template <typename MultiAction, unsigned MaxSteps>
struct multi_action_helper<MultiAction, MaxSteps, 1>
{
    static unsigned const step= MaxSteps - 1;

    void operator() (MultiAction const& action) const
    {
	action(step);
    }

    void operator() (MultiAction& action) const
    {
	action(step);
    }    
};


template <typename MultiAction, unsigned Steps>
struct multi_action_block
{
    void operator() (MultiAction const& action) const
    {
	multi_action_helper<MultiAction, Steps, Steps>()(action);
    }

    void operator() (MultiAction& action) const
    {
	multi_action_helper<MultiAction, Steps, Steps>()(action);
    }
};



} // namespace mtl

#endif // MTL_MULTI_ACTION_BLOCK_INCLUDE
