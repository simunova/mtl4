concept Associative<typename Operation, typename Element>
{
    axiom Associativity(Operation op, Element x, Element y, Element z)
    {
	op(x, op(y, z)) == op(op(x, y), z); 
    }
};
    
concept C1<typename Operation, typename Element>
: Associative<Operation, Element>
{
    axiom A1(Operation op, Element x) { }
};

concept C2<typename Operation, typename Element>
: Associative<Operation, Element>
{
    axiom A2(Operation op, Element x) { }
};

concept ConceptDiamond<typename Operation, typename Element>
: C1<Operation, Element>,
  C2<Operation, Element>
{};


int main()
{
    return 0;
}
