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

#ifndef MTL_CONTIGUOUS_MEMORY_BLOCK_INCLUDE
#define MTL_CONTIGUOUS_MEMORY_BLOCK_INCLUDE

#include <cassert>
#include <algorithm>
#include <boost/static_assert.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/assert.hpp>
#include <boost/numeric/mtl/matrix/dimension.hpp>
#include <boost/numeric/mtl/detail/index.hpp>
#include <boost/numeric/mtl/operation/clone.hpp>



namespace mtl { namespace detail {
using std::size_t;
  
// Macro MTL_ENABLE_ALIGNMENT is by default not set

// Minimal size of memory allocation using alignment
#ifndef MTL_ALIGNMENT_LIMIT
#  define MTL_ALIGNMENT_LIMIT 1024
#endif

// Alignment in memory
#ifndef MTL_ALIGNMENT
#  define MTL_ALIGNMENT 128
#endif


// Size helper for static size
template <unsigned Size>
struct size_helper
{
    typedef size_helper self;

    size_helper() {}
    explicit size_helper(std::size_t size)
    {
	set_size(size);
    }

# ifndef MTL_IGNORE_STATIC_SIZE_VIOLATION
    void set_size(std::size_t MTL_DEBUG_ARG(size))
    {	MTL_CRASH_IF(Size != size, "Try to change a static size!"); }
# else
    void set_size(std::size_t) {}
# endif

    std::size_t used_memory() const { return Size;  }
    friend void swap(self&, self&) {}
};

// Manage size only if template parameter is 0, i.e. dynamic size
template <>
struct size_helper<0>
{
    typedef size_helper self;

    size_helper(std::size_t size= 0) : my_used_memory(size) {}

    void set_size(std::size_t size)
    {
	my_used_memory= size;
    }

    std::size_t used_memory() const
    {
	return my_used_memory;
    }

    friend void swap(self& x, self& y) 
    {
	std::swap(x.my_used_memory, y.my_used_memory);
    }

  protected:
    std::size_t                               my_used_memory;
};


// Encapsulate behavior of alignment

# ifdef MTL_ENABLE_ALIGNMENT

    template <typename Value>
    struct alignment_helper
    {
	typedef alignment_helper self;

	alignment_helper() : malloc_address(0) {}

	Value* alligned_alloc(std::size_t size)
	{
	    if (size == 0)
		return 0;

	    bool        align= size * sizeof(value_type) >= MTL_ALIGNMENT_LIMIT;
	    std::size_t bytes= size * sizeof(value_type);
	    
	    if (align)
		bytes+= MTL_ALIGNMENT - 1;

	    char* p= malloc_address= new char[bytes];
	    if (align)
		while ((long int)(p) % MTL_ALIGNMENT) p++;
		// p+= MTL_ALIGNMENT - (long int)(p) % MTL_ALIGNMENT;

	    return reinterpret_cast<value_type*>(p);
	}

	void aligned_delete(bool is_own, Value*& data)
	{
	    if (is_own && malloc_address) delete[] malloc_address;
	    data= 0;
	}

	friend void swap(self& x, self& y) 
	{
	    using std::swap
	    swap(x.malloc_address, y.malloc_address);
	}

      private:	    
	char*                                     malloc_address;
    };

# else

    template <typename Value>
    struct alignment_helper
    {
	typedef alignment_helper self;

	Value* alligned_alloc(std::size_t size)	{  return size > 0 ? new Value[size] : (Value*)(0); }

	void aligned_delete(bool is_own, Value*& data)
	{
	    if (is_own && data != 0) // std::cout << "Delete " << data << '\n', 
		delete[] data, data= 0;
	}

	friend void swap(self&, self&) {}
    };

# endif


template <typename Value, bool OnStack, unsigned Size> // for data on stack
struct memory_crtp
//    : public contiguous_memory_block<Value, OnStack, Size>
{
    typedef contiguous_memory_block<Value, OnStack, Size> base;

    static bool const                         on_stack= OnStack;
    
    typedef Value                             value_type;
    typedef value_type*                       pointer_type;
    typedef const value_type*                 const_pointer_type;

    // offset of key (pointer) w.r.t. data 
    // values must be stored consecutively
    size_t offset(const Value* p) const 
    { 
      return p - static_cast<const base&>(*this).data; 
    }

    // returns pointer to data
    pointer_type elements()
    {
      return static_cast<base&>(*this).data; 
    }

    // returns const pointer to data
    const_pointer_type elements() const 
    {
      return static_cast<const base&>(*this).data; 
    }

    // returns n-th value in consecutive memory
    // (whatever this means in the corr. matrix format)
    value_type& value_n(size_t offset)
    { 
      return static_cast<base&>(*this).data[offset]; 
    }

    // returns n-th value in consecutive memory
    // (whatever this means in the corr. matrix format)
    const value_type& value_n(size_t offset) const 
    { 
      return static_cast<const base&>(*this).data[offset]; 
    }
    
};

// OnStack == false -> data on heap 
template <typename Value, bool OnStack, unsigned Size>
struct contiguous_memory_block
    : public size_helper<Size>,
      public alignment_helper<Value>,
      public memory_crtp<Value, OnStack, Size>
{
    typedef Value                             value_type;
    typedef contiguous_memory_block           self;
    typedef size_helper<Size>                 size_base;
    typedef alignment_helper<Value>           alignment_base;
    typedef memory_crtp<Value, OnStack, Size> crtp_base;

    /// Category of memory, determines behaviour
    enum c_t {own,         //< My own memory: allocate and free it
	      external,    //< Memory, complete memory block of other item, only reference 
	      view         //< View of other's memory (e.g. sub-matrix), different construction than external
    };

  private:

    void alloc(std::size_t size)
    {
	category= own;
	this->set_size(size);
	data= this->alligned_alloc(this->used_memory());
    }

    void delete_it()
    {
	this->aligned_delete(category == own, data);
    }

    template <typename Other>
    void copy_construction(const Other& other)
    {
	using std::copy;
	category= own;
	// std::cout << "Copied in copy constructor.\n";	
	alloc(other.used_memory());
	// std::cout << "My address: " << data << ", other address: " << other.data << '\n';
	copy(other.data, other.data + other.used_memory(), data);
    }

    void move_construction(self& other)
    {
	// std::cout << "Data moved in constructor.\n";
	category= own; data= 0;
	swap(*this, other);
    }

    // Copy the arguments of a view (shallowly) and leave original as it is
    void copy_view(const self& other)
    {
	// std::cout << "View copied (shallowly).\n";
	assert(other.category == view);
	category= view;
	this->set_size(other.used_memory());
	data= other.data;
    }

    template <typename Other>
    void copy_assignment(const Other& other)
    {
	// std::cout << "Copied in assignment.\n";	
	if (this->used_memory() == 0)
	    alloc(other.used_memory());
	MTL_CRASH_IF(this->used_memory() != other.used_memory(), "Incompatible size");
	std::copy(other.data, other.data + other.used_memory(), data);
    }

  public:
    contiguous_memory_block() : category(own), data(0) {}

    explicit contiguous_memory_block(Value *data, std::size_t size, bool is_view= false) 
	: size_base(size), category(is_view ? view : external), data(data)
    {}    

    explicit contiguous_memory_block(std::size_t size) : category(own)
    {
	// std::cout << "Constructor with size.\n";
	alloc(size);
	// std::cout << "New block at " << data << '\n';
    }

    // Default copy constructor
    contiguous_memory_block(const self& other) : size_base(other)
    {
	// std::cout << "Copy constructor (same type).\n";	
	if (other.category == view)
	    copy_view(other);
	else
	    copy_construction(other);
    }

#ifdef MTL_WITH_MOVE
    contiguous_memory_block(self&& other) // : data(other.data), category(other.category) {}
    {
	move_construction(other);
    }
#endif

    // Force copy construction
    contiguous_memory_block(const self& other, clone_ctor)
    {
	// std::cout << "(Forced) Copy constructor (same type).\n";	
	copy_construction(other);
    }

    // Other types must be copied always
    template<typename Value2, bool OnStack2, unsigned Size2>
    explicit contiguous_memory_block(const contiguous_memory_block<Value2, OnStack2, Size2>& other)
    {
	// std::cout << "Copy constructor (different type).\n";
	copy_construction(other);
    }

#ifdef MTL_WITH_MOVE
    self& operator=(self&& other)
    {
	move_assignment(other);
	return *this;
    }

    self& operator=(const self& other)
    {
	copy_assignment(other);
	return *this;
    }
#elif defined(MTL_MEMORY_BLOCK_MOVE_EMULATION)
    // Operator takes parameter by value and consumes it
    self& operator=(self other)
    {
	move_assignment(other);
	return *this;
    }
#else
    self& operator=(self other)
    {
	copy_assignment(other);
	return *this;
    }
#endif

    // Same behavior as consuming assignment, to be used by derived classes
protected:
    void move_assignment(self& other)
    {
	// std::cout << "Consuming assignment operator (if same type).\n";
	if (category == own && other.category == own)
	    swap(*this, other);
	else
	    copy_assignment(other);
    }

public:
    template<typename Value2, bool OnStack2, unsigned Size2>
    self& operator=(const contiguous_memory_block<Value2, OnStack2, Size2>& other)
    {
	// std::cout << "Assignment from different array type -> Copy.\n";
	copy_assignment(other);
	return *this;
    }


    void set_view() { category= view; }

    void realloc(std::size_t size)
    {
	if (Size == 0) {

	    // If already have memory of the right size we can keep it
	    if (size == this->used_memory()) 
		return;
	    MTL_CRASH_IF(category != own, 
		      "Can't change the size of collections with external memory");
	    delete_it();
	    alloc(size);
	} else {
	    MTL_CRASH_IF(size != Size, "Can't change static size"); 
	}
    }

    ~contiguous_memory_block()
    {
	//std::cout << "Delete block with address " << data << '\n';
	delete_it();
    }

    friend void swap(self& x, self& y)
    {
	using std::swap;
	swap(x.category, y.category);
	std::swap(x.data, y.data);
	swap(static_cast<size_base&>(x), static_cast<size_base&>(y));
	swap(static_cast<alignment_base&>(x), static_cast<alignment_base&>(y));
    }	

protected:
    enum c_t                                  category;
public:
    Value                                     *data;
};

// OnStack == true 
template <typename Value, unsigned Size>
struct contiguous_memory_block<Value, true, Size>
    : public alignment_helper<Value>,
      public memory_crtp<Value, true, Size>
{
    typedef Value                             value_type;
    typedef contiguous_memory_block           self;
    //static bool const                         on_stack= true;

    Value    data[Size];

# ifdef NDEBUG
    contiguous_memory_block() {} // default constructor in release mode
    explicit contiguous_memory_block(std::size_t) {}
# else 
    explicit contiguous_memory_block(std::size_t size= Size)
    {
	MTL_CRASH_IF(Size != size, "Incompatible size!");
    }
# endif

    // Move-semantics ignored for arrays on stack
    contiguous_memory_block(const self& other)
    {
	// std::cout << "Copied in copy constructor (same type).\n";
	std::copy(other.data, other.data+Size, data);
    }


    template<typename Value2, bool OnStack2, unsigned Size2>
    explicit contiguous_memory_block(const contiguous_memory_block<Value2, OnStack2, Size2>& other)
    {
	// std::cout << "Copied in copy constructor (different type).\n";	
	MTL_CRASH_IF(Size != other.used_memory(), "Incompatible size!");
	std::copy(other.data, other.data + other.used_memory(), data);
    }

    self& operator=(const self& other)
    {
	// std::cout << "Assignment from same type.\n";
	std::copy(other.data, other.data+Size, data);
	return *this;
    }

    // For consistency with non-static blocks, to be used by derived classes
protected:
    void move_assignment(self& other)
    {
	std::copy(other.data, other.data+Size, data);
    }

public:
    template<typename Value2, bool OnStack2, unsigned Size2>
    self& operator=(const contiguous_memory_block<Value2, OnStack2, Size2>& other)
    {
	// std::cout << "Assignment from different type.\n";
	MTL_CRASH_IF(Size != other.used_memory(), "Incompatible size!");
	std::copy(other.data, other.data + other.used_memory(), data);
	return *this;
    }


    void realloc(std::size_t MTL_DEBUG_ARG(s)) 
    {
	// Arrays on stack cannot be reallocated but if the size isn't changed we are fine
	assert(s == Size); 
    }

    std::size_t used_memory() const
    {
	return Size;
    }

  protected:
    enum c_t {own};
    static const c_t category= own;
};


}} // namespace mtl::detail

namespace mtl {
    template <typename Value, bool OnStack, unsigned Size>
    struct is_clonable< detail::contiguous_memory_block<Value, OnStack, Size> > : boost::mpl::bool_<!OnStack> {};
}

#endif // MTL_CONTIGUOUS_MEMORY_BLOCK_INCLUDE
