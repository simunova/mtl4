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

#ifndef MTL_VECTOR_MAP_VIEW_INCLUDE
#define MTL_VECTOR_MAP_VIEW_INCLUDE

#include <utility>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/utility/copy_expression_const_ref_container.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>
#include <boost/numeric/mtl/operation/tfunctor.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/mtl/operation/real.hpp>
#include <boost/numeric/mtl/operation/imag.hpp>
#include <boost/numeric/mtl/operation/signum.hpp>
#include <boost/numeric/mtl/vector/vec_expr.hpp>


namespace mtl { namespace vec { namespace detail {
    // Forward declaration for friend declaration
    template <typename, typename> struct map_value;
}}}

namespace mtl { namespace vec {

template <typename Functor, typename Vector> 
struct map_view 
  : public vec_expr< map_view<Functor, Vector> >
{
    typedef map_view                                   self;
    typedef vec_expr< self >                           expr_base;
    typedef Vector                                     other;

    typedef typename Functor::result_type              value_type;
    typedef typename Functor::result_type              const_reference;
    typedef typename Collection<Vector>::size_type     size_type;

    map_view (const Functor& functor, const other& ref) 
      : expr_base(*this), functor(functor), ref(ref) 
    {
  ref.delay_assign();
    }
    
    map_view (const Functor& functor, boost::shared_ptr<Vector> p) 
      : expr_base(*this), functor(functor), my_copy(p), ref(*p)
    {
  ref.delay_assign();
    }

#ifdef MTL_WITH_MOVE    
  map_view (self&& that) : my_copy(std::move(that.my_copy)), functor(that.functor), ref(that.ref) {}
  map_view (const self& that) : functor(that.functor), ref(that.ref) { assert(that.my_copy.use_count() == 0); }
#endif

    friend size_type inline num_rows(const self& v) { return num_rows(v.ref); }
    friend size_type inline num_cols(const self& v) { return num_cols(v.ref); }

    size_type stride() const {   return ref.stride(); }
    const_reference operator() (size_type i) const { return functor(ref(i)); }
    const_reference operator[] (size_type i) const { return functor(ref[i]); }
    void delay_assign() const {}
    
    template <typename, typename> friend struct detail::map_value;

  protected:
    boost::shared_ptr<Vector>           my_copy;
  public:
    Functor           functor;
    // ref is a const& if Vector is a true vector and a copy if it is an expression
    typename mtl::traits::copy_expression_const_ref_container<Vector>::type ref;
};


template <typename Functor, typename Vector> 
inline std::size_t size(const map_view<Functor, Vector>& v)
{     return mtl::vec::size(v.ref); }

// ================
// Free functions
// ================


    namespace detail {

  template <typename Functor, typename Vector> 
  struct map_value
  {
      typedef typename Vector::key_type                      key_type;
      typedef typename vec::map_view<Functor, Vector>::value_type value_type;
      
      map_value(vec::map_view<Functor, Vector> const& map_vector) 
    : map_vector(map_vector), its_value(map_vector.ref) 
      {}

      value_type operator() (key_type const& key) const
      {
    return map_vector.functor(its_value(key));
      }

    protected:
      vec::map_view<Functor, Vector> const&   map_vector;
      typename ::mtl::traits::const_value<Vector>::type its_value;
        };

    } // detail

}} // namespace mtl::vector



namespace mtl { namespace traits {

    // ================
    // Property maps
    // ================

    template <typename Functor, typename Vector> 
    struct index<vec::map_view<Functor, Vector> >
  : public index<Vector>
    {};

    template <typename Functor, typename Vector> 
    struct const_value<vec::map_view<Functor, Vector> >
    {
  typedef vec::detail::map_value<Functor, Vector>  type;
    };


    // ================
    // Range generators
    // ================

    // Use range_generator of original vector
    template <typename Tag, typename Functor, typename Vector> 
    struct range_generator<Tag, vec::map_view<Functor, Vector> >
  : public range_generator<Tag, Vector>
    {};

}} // mtl::traits

namespace mtl { namespace vec {

template <typename Scaling, typename Vector>
struct scaled_view
    : public map_view<tfunctor::scale<Scaling, typename Vector::value_type>, Vector>
{
    typedef tfunctor::scale<Scaling, typename Vector::value_type>  functor_type;
    typedef map_view<functor_type, Vector>                         base;
    typedef scaled_view                                            self;

    explicit scaled_view(const Scaling& scaling, const Vector& vector)
      : base(functor_type(scaling), vector)
    {}
    
    explicit scaled_view(const Scaling& scaling, boost::shared_ptr<Vector> p)
      : base(functor_type(scaling), p)
    {}

#ifdef MTL_WITH_MOVE    
    scaled_view (self&& that) : base(that) {}
    scaled_view (const self& that) : base(that) {}
#endif
};

// added by Hui Li
template <typename Vector, typename RScaling>
struct rscaled_view
  : public map_view<tfunctor::rscale<typename Vector::value_type, RScaling>, Vector>
{
    typedef tfunctor::rscale<typename Vector::value_type, RScaling>  functor_type;
    typedef map_view<functor_type, Vector>                          base;
    typedef rscaled_view                                            self;
  
    explicit rscaled_view(const Vector& vector, const RScaling& rscaling)
      : base(functor_type(rscaling), vector)
    {}
  
    explicit rscaled_view(boost::shared_ptr<Vector> p, const RScaling& rscaling)
      : base(functor_type(rscaling), p)
    {}

#ifdef MTL_WITH_MOVE    
    rscaled_view (self&& that) : base(that) {}
    rscaled_view (const self& that) : base(that) {}
#endif
};
  

// added by Hui Li
template <typename Vector, typename Divisor>
struct divide_by_view
  : public map_view<tfunctor::divide_by<typename Vector::value_type, Divisor>, Vector>
{
    typedef tfunctor::divide_by<typename Vector::value_type, Divisor>  functor_type;
    typedef map_view<functor_type, Vector>                             base;
    typedef divide_by_view                                             self;
  
    explicit divide_by_view(const Vector& vector, const Divisor& div)
      : base(functor_type(div), vector)
    {}
  
    explicit divide_by_view(boost::shared_ptr<Vector> p, const Divisor& div)
      : base(functor_type(div), p)
    {}
  
#ifdef MTL_WITH_MOVE    
    divide_by_view (self&& that) : base(that) {}
    divide_by_view (const self& that) : base(that) {}
#endif
};
  
/// View for raising vector element to power of scalar exponent
template <typename Vector, typename Exponent>
struct pow_by_view
  : public map_view<tfunctor::pow_by<typename Vector::value_type, Exponent>, Vector>
{
    typedef tfunctor::pow_by<typename Vector::value_type, Exponent>  functor_type;
    typedef map_view<functor_type, Vector>                           base;
    typedef pow_by_view                                              self;
        
    explicit pow_by_view(const Vector& vector, const Exponent& exp)
      : base(functor_type(exp), vector)
    {}
        
    explicit pow_by_view(boost::shared_ptr<Vector> p, const Exponent& exp)
      : base(functor_type(exp), p)
    {}
        
#ifdef MTL_WITH_MOVE    
    pow_by_view (self&& that) : base(that) {}
    pow_by_view (const self& that) : base(that) {}
#endif
};


template <typename Vector>
struct conj_view
  : public map_view<mtl::sfunctor::conj<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::conj<typename Vector::value_type>            functor_type;
    typedef map_view<functor_type, Vector>                         base;
    typedef conj_view                                                   self;

    explicit conj_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit conj_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    conj_view (self&& that) : base(that) {}
    conj_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct real_view
  : public map_view<mtl::sfunctor::real<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::real<typename Vector::value_type>            functor_type;
    typedef map_view<functor_type, Vector>                         base;
  typedef real_view                                              self;

    explicit real_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit real_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    real_view (self&& that) : base(that) {}
    real_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct imag_view
  : public map_view<mtl::sfunctor::imag<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::imag<typename Vector::value_type>            functor_type;
    typedef map_view<functor_type, Vector>                         base;
    typedef imag_view                                              self;

    explicit imag_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit imag_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    imag_view (self&& that) : base(that) {}
    imag_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct negate_view
  : public map_view<mtl::sfunctor::negate<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::negate<typename Vector::value_type>            functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef negate_view                                                   self;

    explicit negate_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit negate_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    negate_view (self&& that) : base(that) {}
    negate_view (const self& that) : base(that) {}
#endif
};


template <typename Vector>
struct abs_view
  : public map_view<mtl::sfunctor::abs<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::abs<typename Vector::value_type>               functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef abs_view                                                      self;

    explicit abs_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit abs_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    abs_view (self&& that) : base(that) {}
    abs_view (const self& that) : base(that) {}
#endif
};

// Inverse trigonometric

template <typename Vector>
struct acos_view
  : public map_view<mtl::sfunctor::acos<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::acos<typename Vector::value_type>              functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef acos_view                                                     self;

    explicit acos_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit acos_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    acos_view (self&& that) : base(that) {}
    acos_view (const self& that) : base(that) {}
#endif
};


# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Vector>
struct acosh_view
  : public map_view<mtl::sfunctor::acosh<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::acosh<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef acosh_view                                                    self;

    explicit acosh_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit acosh_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    acosh_view (self&& that) : base(that) {}
    acosh_view (const self& that) : base(that) {}
#endif
};
#endif

template <typename Vector>
struct asin_view
  : public map_view<mtl::sfunctor::asin<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::asin<typename Vector::value_type>              functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef asin_view                                                     self;

    explicit asin_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit asin_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    asin_view (self&& that) : base(that) {}
    asin_view (const self& that) : base(that) {}
#endif
};


# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Vector>
struct asinh_view
  : public map_view<mtl::sfunctor::asinh<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::asinh<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef asinh_view                                                    self;

    explicit asinh_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit asinh_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    asinh_view (self&& that) : base(that) {}
    asinh_view (const self& that) : base(that) {}
#endif
};
#endif

template <typename Vector>
struct atan_view
  : public map_view<mtl::sfunctor::atan<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::atan<typename Vector::value_type>              functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef atan_view                                                     self;

    explicit atan_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit atan_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    atan_view (self&& that) : base(that) {}
    atan_view (const self& that) : base(that) {}
#endif
};


# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Vector>
struct atanh_view
  : public map_view<mtl::sfunctor::atanh<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::atanh<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef atanh_view                                                    self;

    explicit atanh_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit atanh_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    atanh_view (self&& that) : base(that) {}
    atanh_view (const self& that) : base(that) {}
#endif
};
#endif


// Trigonometric

template <typename Vector>
struct cos_view
  : public map_view<mtl::sfunctor::cos<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::cos<typename Vector::value_type>              functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef cos_view                                                     self;

    explicit cos_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit cos_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    cos_view (self&& that) : base(that) {}
    cos_view (const self& that) : base(that) {}
#endif
};


template <typename Vector>
struct cosh_view
  : public map_view<mtl::sfunctor::cosh<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::cosh<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef cosh_view                                                    self;

    explicit cosh_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit cosh_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    cosh_view (self&& that) : base(that) {}
    cosh_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct sin_view
  : public map_view<mtl::sfunctor::sin<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::sin<typename Vector::value_type>              functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef sin_view                                                     self;

    explicit sin_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit sin_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    sin_view (self&& that) : base(that) {}
    sin_view (const self& that) : base(that) {}
#endif
};


template <typename Vector>
struct sinh_view
  : public map_view<mtl::sfunctor::sinh<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::sinh<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef sinh_view                                                    self;

    explicit sinh_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit sinh_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    sinh_view (self&& that) : base(that) {}
    sinh_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct tan_view
  : public map_view<mtl::sfunctor::tan<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::tan<typename Vector::value_type>              functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef tan_view                                                     self;

    explicit tan_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit tan_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    tan_view (self&& that) : base(that) {}
    tan_view (const self& that) : base(that) {}
#endif
};


template <typename Vector>
struct tanh_view
  : public map_view<mtl::sfunctor::tanh<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::tanh<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef tanh_view                                                    self;

    explicit tanh_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit tanh_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    tanh_view (self&& that) : base(that) {}
    tanh_view (const self& that) : base(that) {}
#endif
};

// Rounding
template <typename Vector>
struct ceil_view
  : public map_view<mtl::sfunctor::ceil<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::ceil<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef ceil_view                                                    self;

    explicit ceil_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit ceil_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    ceil_view (self&& that) : base(that) {}
    ceil_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct floor_view
  : public map_view<mtl::sfunctor::floor<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::floor<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef floor_view                                                    self;

    explicit floor_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit floor_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    floor_view (self&& that) : base(that) {}
    floor_view (const self& that) : base(that) {}
#endif
};

# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Vector>
struct round_view
  : public map_view<mtl::sfunctor::round<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::round<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef round_view                                                    self;

    explicit round_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit round_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    round_view (self&& that) : base(that) {}
    round_view (const self& that) : base(that) {}
#endif
};
#endif

# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Vector>
struct trunc_view
  : public map_view<mtl::sfunctor::trunc<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::trunc<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef trunc_view                                                    self;

    explicit trunc_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit trunc_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    trunc_view (self&& that) : base(that) {}
    trunc_view (const self& that) : base(that) {}
#endif
};
#endif

template <typename Vector>
struct log_view
  : public map_view<mtl::sfunctor::log<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::log<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef log_view                                                    self;

    explicit log_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit log_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    log_view (self&& that) : base(that) {}
    log_view (const self& that) : base(that) {}
#endif
};

# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Vector>
struct log2_view
  : public map_view<mtl::sfunctor::log2<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::log2<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef log2_view                                                    self;

    explicit log2_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit log2_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    log2_view (self&& that) : base(that) {}
    log2_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct log10_view
  : public map_view<mtl::sfunctor::log10<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::log10<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef log10_view                                                    self;

    explicit log10_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit log10_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    log10_view (self&& that) : base(that) {}
    log10_view (const self& that) : base(that) {}
#endif
};
# endif

template <typename Vector>
struct exp_view
  : public map_view<mtl::sfunctor::exp<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::exp<typename Vector::value_type>               functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef exp_view                                                      self;

    explicit exp_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit exp_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    exp_view (self&& that) : base(that) {}
    exp_view (const self& that) : base(that) {}
#endif
};

# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Vector>
struct exp2_view
  : public map_view<mtl::sfunctor::exp2<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::exp2<typename Vector::value_type>              functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef exp2_view                                                     self;

    explicit exp2_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit exp2_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    exp2_view (self&& that) : base(that) {}
    exp2_view (const self& that) : base(that) {}
#endif
};
#endif

template <typename Vector>
struct exp10_view
  : public map_view<mtl::sfunctor::exp10<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::exp10<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef exp10_view                                                    self;

    explicit exp10_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit exp10_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    exp10_view (self&& that) : base(that) {}
    exp10_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct sqrt_view
  : public map_view<mtl::sfunctor::sqrt<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::sqrt<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef sqrt_view                                                    self;

    explicit sqrt_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit sqrt_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    sqrt_view (self&& that) : base(that) {}
    sqrt_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct rsqrt_view
  : public map_view<mtl::sfunctor::rsqrt<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::rsqrt<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef rsqrt_view                                                    self;

    explicit rsqrt_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit rsqrt_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    rsqrt_view (self&& that) : base(that) {}
    rsqrt_view (const self& that) : base(that) {}
#endif
};

# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Vector>
struct erf_view
  : public map_view<mtl::sfunctor::erf<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::erf<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef erf_view                                                    self;

    explicit erf_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit erf_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    erf_view (self&& that) : base(that) {}
    erf_view (const self& that) : base(that) {}
#endif
};

template <typename Vector>
struct erfc_view
  : public map_view<mtl::sfunctor::erfc<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::erfc<typename Vector::value_type>             functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef erfc_view                                                    self;

    explicit erfc_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit erfc_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    erfc_view (self&& that) : base(that) {}
    erfc_view (const self& that) : base(that) {}
#endif
};
# endif

template <typename Vector>
struct signum_view
  : public map_view<mtl::sfunctor::signum<typename Vector::value_type>, Vector>
{
    typedef mtl::sfunctor::signum<typename Vector::value_type>              functor_type;
    typedef map_view<functor_type, Vector>                                base;
    typedef signum_view                                                     self;

    explicit signum_view(const Vector& vector)
      : base(functor_type(), vector)
    {}
    
    explicit signum_view(boost::shared_ptr<Vector> p)
      : base(functor_type(), p)
    {}

#ifdef MTL_WITH_MOVE    
    signum_view (self&& that) : base(that) {}
    signum_view (const self& that) : base(that) {}
#endif
};

}} // namespace mtl::vector

namespace mtl { namespace sfunctor {


    template <typename Vector>
    struct conj_aux<Vector, tag::vector>
    {
  typedef mtl::vec::conj_view<Vector> result_type;

  static inline result_type apply(const Vector& vector)
  {
      return result_type(vector);
  }

  result_type operator() (const Vector& vector) const
  {
      return apply(vector);
  }
    };

}}




#endif // MTL_VECTOR_MAP_VIEW_INCLUDE
