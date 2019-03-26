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

#ifndef MTL_POSITIVE_REAL_INCLUDE
#define MTL_POSITIVE_REAL_INCLUDE

namespace mtl {

class positive_real 
{
protected:
    double value; 
public:
    positive_real(double m): value(m) 
    {
	if (m < 0.0) throw "Negative value not allowed!\n"; 
    }

    double get_value() const 
    {
	return value; 
    }

    positive_real operator+(positive_real const& y) const
    {
	return value + y.value;
    }

    positive_real& operator+=(positive_real const& y)
    {
	value+= y.value;
	return *this;
    }

    positive_real operator*(positive_real const& y) const
    {
	return value * y.value;
    }

    positive_real& operator*=(positive_real const& y)
    {
	value*= y.value;
	return *this;
    }

    positive_real operator/(positive_real const& y) const
    {
	return value / y.value;
    }

    positive_real& operator/=(positive_real const& y)
    {
	value/= y.value;
	return *this;
    }

    bool operator==(positive_real const& y) const
    {
	return value == y.value;
    }
    
    bool operator!=(positive_real const& y) const
    {
	return value != y.value;
    }
};
 
 
inline std::ostream& operator<< (std::ostream& stream, const positive_real& a) 
{
    return stream << a.get_value(); 
}

struct checked_positive_real : public positive_real 
{
    void check_range(checked_positive_real const& x)
    {
#     ifndef NDEBUG
	if (x.get_value() == 0) throw "In struct group_real: x is zero";
	if (std::isinf(x.get_value()))  throw "In struct group_real: x is infinity";
#     endif
    }

    void check_range() 
    {
	check_range(*this);
    }

    explicit checked_positive_real(float x) : positive_real(x)
    {
	check_range();
    }

    checked_positive_real operator*(checked_positive_real const& y) const
    {
	checked_positive_real tmp(*this);
	tmp*= y;
	return tmp;
    }

    checked_positive_real& operator*=(checked_positive_real const& y)
    {
	value*= y.get_value();
	check_range();
	return *this;
    }

    checked_positive_real operator/(checked_positive_real const& y) const
    {
	checked_positive_real tmp(*this);
	tmp/= y;
	return tmp;
    }

    checked_positive_real& operator/=(checked_positive_real const& y)
    {
	value/= y.get_value();
	check_range();
	return *this;
    }
    
};


} // namespace mtl



#endif // MTL_POSITIVE_REAL_INCLUDE
