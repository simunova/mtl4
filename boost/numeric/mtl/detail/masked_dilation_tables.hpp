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

#ifndef MTL_MASKED_DILATION_TABLES_INCLUDE
#define MTL_MASKED_DILATION_TABLES_INCLUDE

#include <iostream>
#include <iomanip>
#include <cassert>

namespace mtl { namespace dilated {

template <class T, T Mask>
struct masked_dilation_tables
{
    typedef masked_dilation_tables     self;
    typedef T                          value_type;
    const static T                     mask= Mask;

    static const T n_bytes= sizeof(T);           // number of bytes of the unmasked value of this type
    typedef T                          lookup_type[n_bytes][256];
    typedef T                          mp_type[n_bytes];
    typedef T                          it_type[n_bytes];  // why int ??? switch to T

  protected:
    static lookup_type*                my_mask_lut, *my_unmask_lut;
    static mp_type*                    my_mask_piece;
    static it_type*                    my_mask_size, *my_mask_shift_table, *my_unmask_shift_table;
    static int                         n_valid_table, instances;
  public:
    
    masked_dilation_tables()
    {
	if (instances++ == 0)
	    compute_tables();
    }

    ~masked_dilation_tables()
    {
	if (--instances == 0) {
	    delete[] my_mask_lut;         
	    delete[] my_unmask_lut;       
	    delete[] my_mask_piece;        
	    delete[] my_mask_size;         
	    delete[] my_mask_shift_table;  
	    delete[] my_unmask_shift_table;
	}
    }

    lookup_type& mask_lut()   { return *my_mask_lut; }
    lookup_type& unmask_lut() { return *my_unmask_lut; }
    mp_type& mask_piece()     { return *my_mask_piece; }
    it_type& mask_size()      {	return *my_mask_size;  }
    it_type& mask_shift_table() { return *my_mask_shift_table; }
    it_type& unmask_shift_table() { return *my_unmask_shift_table; }

private:

    // get mask of the style 0xfff...
    static T get_f_mask(T n_bits) { return (1 << n_bits) - 1;   }

    T inc(T i, T mask) { return ((i - mask) & mask);    }

    void compute_tables() 
    {
	// std::cout << "computing tables! " << std::endl;
	init();

	// compute the mask table
	for (int j = 0; j < n_valid_table; ++j) {
	    T f_mask = get_f_mask(mask_size()[j]), i, ii; 
	    for (i = 0, ii = 0; i < 256; ++i, ii = inc(ii, mask_piece()[j])) 
		mask_lut()[j][i] =  (ii & f_mask) << mask_shift_table()[j]; // need to shift 
	}

	// compute the unmask table
	T f_mask = get_f_mask(8);
	for (T j = 0; j < sizeof(T); ++j) {
	    T t_mask = (Mask >> (8*j)) & f_mask, i, ii;
	    for(i = 0, ii = 0; ii < t_mask; ii = inc(ii, t_mask), ++i) 
		unmask_lut()[j][ii] =  i << unmask_shift_table()[j];
	    // set the value for the last one
	    unmask_lut()[j][t_mask] =  i << unmask_shift_table()[j];       
	}
    }

    void allocate()
    {
	my_mask_lut=            new lookup_type[1];
	my_unmask_lut=          new lookup_type[1];
	my_mask_piece=          new mp_type[1];
	my_mask_size=           new it_type[1];
	my_mask_shift_table=    new it_type[1];
	my_unmask_shift_table=  new it_type[1];
    }


    // initialize needed parameters
    void init() 
    {
	allocate();
	assert(count_n_ones(Mask) > 0);
	n_valid_table= (count_n_ones(Mask) + 7) / 8; // calculate the number of valid table
	set_mask();    
    }

    // return the number of 1's in the mask
    int count_n_ones(T t) 
    {
	int n_ones = 0;
	for (; t; t>>= 1)
	    if(t & 1) ++n_ones;
	return n_ones;
    }

    // return the number of valid bits in the mask
    int count_bits(T t) 
    {
	int bits = 0;
	for (; t; t>>= 1)
	    ++bits;
	return bits;
    }


    // set mask pieces
    void set_mask() 
    {
	// set the unmask shift table
	unmask_shift_table()[0] = 0;
	T t_mask = Mask, tmp, count;
	for (T i = 1; i < n_bytes; ++i) {
	    tmp = t_mask & get_f_mask(8);
	    count = count_n_ones(tmp);
	    unmask_shift_table()[i] = count + unmask_shift_table()[i - 1];
	    t_mask >>= 8;
	}

	mask_shift_table()[0] = 0;  // don't need shift for the first table
	// if there is only 8 or less 1's in the mask,
	// only one table is needed
	if (n_valid_table == 1) {
	    mask_piece()[0] = Mask; 
	    mask_size()[0] = count_bits(Mask);
	    return;
	}

	t_mask = Mask;
	for (int i = 0; i < n_valid_table - 1; ++i) {
	    T n_bits = 0, tmp = t_mask;
	    for (T n_ones= 0; n_ones < 8; ++n_bits) {
		if ((t_mask & 0x01) == 1) ++n_ones;
		t_mask = t_mask >>1;
	    } 
	    // set the ith piece of mask, which must contains 8 1's
	    mask_piece()[i] = get_f_mask(n_bits) & tmp;
	    
	    // set the mask size table
	    mask_size()[i] = n_bits;
	    
	    // set shift table
	    mask_shift_table()[i + 1] = n_bits + mask_shift_table()[i];
	}

	// set the last piece of mask, which may contain less than 8 1's
	// set the number of bits of the last mask
	mask_piece()[n_valid_table - 1 ] = t_mask;    
	mask_size()[n_valid_table - 1] = count_bits(t_mask);
    }

public:
    
    // convert to masked integer
    T to_masked(T x) 
    {
	T result = 0;
	for (int i = 0; i < n_valid_table; ++i)
	    result += mask_lut()[i][0xff & (x >> (8*i)) ];
	return result;
    }


    // convert to unmasked integer
    T to_unmasked(T x) 
    {
	T result = 0;
	x &= Mask;
	for (T i = 0; i < n_bytes; ++i) {
	    result += unmask_lut()[i][0xff & (x >> (8*i)) ];
	}
	return result;
    }
};

template <class T, T Mask>
typename masked_dilation_tables<T, Mask>::lookup_type* masked_dilation_tables<T, Mask>::my_mask_lut= 0;

template <class T, T Mask>
typename masked_dilation_tables<T, Mask>::lookup_type* masked_dilation_tables<T, Mask>::my_unmask_lut= 0;

template <class T, T Mask>
typename masked_dilation_tables<T, Mask>::mp_type* masked_dilation_tables<T, Mask>::my_mask_piece= 0;

template <class T, T Mask>
typename masked_dilation_tables<T, Mask>::it_type* masked_dilation_tables<T, Mask>::my_mask_size= 0;

template <class T, T Mask>
typename masked_dilation_tables<T, Mask>::it_type* masked_dilation_tables<T, Mask>::my_mask_shift_table= 0;

template <class T, T Mask>
typename masked_dilation_tables<T, Mask>::it_type* masked_dilation_tables<T, Mask>::my_unmask_shift_table= 0;

template <class T, T Mask>
int masked_dilation_tables<T, Mask>::n_valid_table= 0;

template <class T, T Mask>
int masked_dilation_tables<T, Mask>::instances= 0;


// Masking: syntax e.g. mask<0x55555555>(7);
// Mask must be in front of T -> need casting :-(
template <std::size_t Mask, typename T>
inline T mask(T const& value)
{
    static masked_dilation_tables<T, T(Mask)>  tables;
    return tables.to_masked(value);
}


// Masking: syntax e.g. mask(7, table_object);
template <typename T, T Mask>
inline T mask(T const& value, masked_dilation_tables<T, Mask> tables)
{
    return tables.to_masked(value);
}


// Unmasking: syntax e.g. unmask<0x55555555>(7);
// Mask must be in front of T -> need casting :-(
template <std::size_t Mask, typename T>
inline T unmask(T const& value)
{
    static masked_dilation_tables<T, T(Mask)>  tables;
    return tables.to_unmasked(value);
}


// Unmasking: syntax e.g. unmask(7, table_object);
template <typename T, T Mask>
inline T unmask(T const& value, masked_dilation_tables<T, Mask> tables)
{
    return tables.to_unmasked(value);
}


// Conversion from Mask1 to Mask2
// syntax e.g. from Morton to Doppler convert<0x55555555, 0x5555ff00>(7); 
// Mask must be in front of T -> need casting :-(
template <long unsigned Mask1, long unsigned Mask2, typename T>
inline T convert(T const& value)
{
    return mask<Mask2>(unmask<Mask1>(value));
}

// Conversion from Mask1 to Mask2
template <long unsigned Mask1, long unsigned Mask2, typename T>
inline T convert(T const& value, masked_dilation_tables<T, Mask1> const& tables1, 
		 masked_dilation_tables<T, Mask2> const& tables2)
{
    return tables2.to_masked(tables1.to_unmasked(value));
}



} // namespace mtl::dilated

  //using dilated::dilated_int;

} // namespace mtl

#endif // MTL_MASKED_DILATION_TABLES_INCLUDE
