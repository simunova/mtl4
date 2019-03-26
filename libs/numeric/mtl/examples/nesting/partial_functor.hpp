template <typename ValueA, typename ParaA, typename ValueB, typename ParaB, 
	  typename ValueC, typename ParaC>
struct mult_ft<dense2D<ValueA, ParaA>, dense2D<ValueB, ParaB>, dense2D<ValueC, ParaC> >
{
    void operator()(dense2D<ValueA, ParaA> const& a, dense2D<ValueB, ParaB> const& b, 
		    dense2D<ValueC, ParaC>& c) 
    { /* Faster code for this set of type triplets ... */ }
};
