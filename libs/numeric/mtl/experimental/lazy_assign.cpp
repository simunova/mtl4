// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschrÃ¤nkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp> 
#include <boost/static_assert.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>

using namespace std; 

template <typename T, typename U, typename Assign>
struct lazy_assign
{
    typedef Assign  assign_type;

    lazy_assign(T& first, const U& second) : first(first), second(second) {} 

    T&       first;
    const U& second;

};

template <typename T, typename U, typename Assign>
void inline evaluate_lazy(lazy_assign<T, U, Assign>& lazy)
{
    Assign::first_update(lazy.first, lazy.second);
}


template <typename T>
struct is_lazy : boost::mpl::false_ {};

template <typename T, typename U, typename Assign>
struct is_lazy<lazy_assign<T, U, Assign> > : boost::mpl::true_ {};


template <typename T>
struct lazy_t
{
    lazy_t(T& data) : data(data) {}

    template <typename U>
    lazy_assign<T, U, mtl::assign::assign_sum> operator=(const U& other) 
    { return lazy_assign<T, U, mtl::assign::assign_sum>(data, other); }

    template <typename U>
    lazy_assign<T, U, mtl::assign::plus_sum> operator+=(const U& other) 
    { return lazy_assign<T, U, mtl::assign::plus_sum>(data, other); }

    template <typename U>
    lazy_assign<T, U, mtl::assign::minus_sum> operator-=(const U& other) 
    { return lazy_assign<T, U, mtl::assign::minus_sum>(data, other); }

    T& data;
};

template <typename T>
inline lazy_t<T> lazy(T& x) 
{ return lazy_t<T>(x); }

template <typename T>
inline lazy_t<const T> lazy(const T& x) 
{ return lazy_t<const T>(x); }

template <typename T>
struct is_vector_reduction : boost::mpl::false_ {};

template <unsigned long Unroll, typename Vector1, typename Vector2, typename ConjOpt>
struct is_vector_reduction<mtl::dot_class<Unroll, Vector1, Vector2, ConjOpt> >
  : boost::mpl::true_ {};

#if 0
template <unsigned long Unroll, typename Vector>
struct is_vector_reduction<mtl::unary_dot_class<Unroll, Vector> >
  : boost::mpl::true_ {};
#endif

template<typename Vector, typename Functor>
struct is_vector_reduction<mtl::lazy_reduction<Vector, Functor> >
  : boost::mpl::true_ {};

template <typename T>
struct index_evaluatable : boost::mpl::false_ {};

template <typename T, typename U, typename Assign>
struct index_evaluatable<lazy_assign<T, U, Assign> >
  : boost::mpl::or_<
      boost::mpl::and_<mtl::traits::is_vector<T>, mtl::traits::is_scalar<U> >,
      boost::mpl::and_<mtl::traits::is_vector<T>, mtl::traits::is_vector<U> >,
      boost::mpl::and_<mtl::traits::is_scalar<T>, is_vector_reduction<U> >
    >
{};

template <typename V1, typename Matrix, typename V2, typename Assign>
struct index_evaluatable<lazy_assign<V1, mtl::mat_cvec_times_expr<Matrix, V2>, Assign> >
  : mtl::traits::is_row_major<Matrix> {};

// Strided traversal would be more expensive than saving from mixed reduction
//  : boost::mpl::or_<mtl::traits::is_row_major<Matrix>, mtl::traits::is_dense<Matrix> > {}; 

template <typename T> struct evaluator_type {};

template <typename T, typename U, typename Assign>
struct evaluator_type<lazy_assign<T, U, Assign> >
  : boost::lazy_enable_if<mtl::traits::is_vector<T>,
			  boost::mpl::if_<mtl::traits::is_vector<U>,
					  mtl::vec_vec_aop_expr<T, U, Assign>, 
					  mtl::vec_scal_aop_expr<T, U, Assign>
					  >
			  >
{};

template <typename T, typename U, typename Assign>
typename boost::enable_if<boost::mpl::and_<mtl::traits::is_vector<T>, mtl::traits::is_vector<U> >, 
			  mtl::vec_vec_aop_expr<T, U, Assign> >::type
inline index_evaluator(lazy_assign<T, U, Assign>& lazy)
{
    return mtl::vec_vec_aop_expr<T, U, Assign>(lazy.first, lazy.second, true);
}


template <typename T, typename U, typename Assign>
typename boost::enable_if<boost::mpl::and_<mtl::traits::is_vector<T>, mtl::traits::is_scalar<U> >, 
			  mtl::vec_scal_aop_expr<T, U, Assign> >::type
inline index_evaluator(lazy_assign<T, U, Assign>& lazy)
{
    return mtl::vec_scal_aop_expr<T, U, Assign>(lazy.first, lazy.second, true);
}

#if 0
template <typename Scalar, typename Vector, typename Assign>
struct unary_dot_index_evaluator
{
    unary_dot_index_evaluator(Scalar& scalar, const Vector& v) 
      : scalar(scalar), v(v) 
    { 
	tmp[0]= tmp[1]= tmp[2]= tmp[3]= Scalar(0); 
    }

    ~unary_dot_index_evaluator() 
    { 
	Scalar s(tmp[0] + tmp[1] + tmp[2] + tmp[3]);
	Assign::apply(scalar, s); 
    }
    
    void operator() (std::size_t i) { mtl::two_norm_functor::update(tmp[0], v[i]); }
    void operator[] (std::size_t i) { (*this)(i); }
    
    template <unsigned Offset>
    void at(std::size_t i) 
    { mtl::two_norm_functor::update(tmp[Offset], v[i+Offset]); }

    Scalar&        scalar;
    Scalar         tmp[4];
    const Vector&  v;
};

template <typename Scalar, typename Vector, typename Assign>
inline std::size_t size(const unary_dot_index_evaluator<Scalar, Vector, Assign>& eval)
{ return size(eval.v); }
#endif


template <typename Scalar, typename Vector, typename Functor, typename Assign>
struct reduction_index_evaluator
{
    reduction_index_evaluator(Scalar& scalar, const Vector& v) 
      : scalar(scalar), v(v) 
    {
	Functor::init(tmp[0]);
	tmp[1]= tmp[2]= tmp[3]= tmp[0];
    }

    ~reduction_index_evaluator() 
    { 
	Functor::finish(tmp[0], tmp[1]);
	Functor::finish(tmp[2], tmp[3]);
	Functor::finish(tmp[0], tmp[2]);
	Assign::apply(scalar, Functor::post_reduction(tmp[0])); // compute sqrt or such if necessary
    }

    template <unsigned Offset>
    void at(std::size_t i) 
    { 
	Functor::update(tmp[Offset], v[i+Offset]); 
    }

    void operator[] (std::size_t i) { at<0>(i); }
    void operator() (std::size_t i) { at<0>(i); }    

    Scalar&        scalar;
    Scalar         tmp[4];
    const Vector&  v;
};

template <typename Scalar, typename Vector, typename Functor, typename Assign>
inline std::size_t size(const reduction_index_evaluator<Scalar, Vector, Functor, Assign>& eval)
{ return size(eval.v); }



template <typename Scalar, typename Vector, typename Functor, typename Assign>
reduction_index_evaluator<Scalar, Vector, Functor, Assign>
inline index_evaluator(lazy_assign<Scalar, mtl::lazy_reduction<Vector, Functor>, Assign>& lazy)
{
    return reduction_index_evaluator<Scalar, Vector, Functor, Assign>(lazy.first, lazy.second.v);
}


template <typename Scalar, typename Vector1, typename Vector2, typename ConjOpt, typename Assign>
struct dot_index_evaluator
{
    dot_index_evaluator(Scalar& scalar, const Vector1& v1, const Vector2& v2) 
      : scalar(scalar), v1(v1), v2(v2) 
    { 
	tmp[0]= tmp[1]= tmp[2]= tmp[3]= Scalar(0); 
    }

    ~dot_index_evaluator() 
    { 
	Scalar s(tmp[0] + tmp[1] + tmp[2] + tmp[3]);
	Assign::apply(scalar, s); 
    }
    
    void operator() (std::size_t i) { tmp[0]+= ConjOpt()(v1[i]) * v2[i]; }
    void operator[] (std::size_t i) { (*this)(i); }

    template <unsigned Offset>
    void at(std::size_t i) 
    { tmp[Offset]+= ConjOpt()(v1[i+Offset]) * v2[i+Offset]; }

    Scalar&        scalar;
    Scalar         tmp[4];
    const Vector1& v1;
    const Vector2& v2;
};

template <typename Scalar, typename Vector1, typename Vector2, typename ConjOpt, typename Assign>
inline std::size_t size(const dot_index_evaluator<Scalar, Vector1, Vector2, ConjOpt, Assign>& eval)
{ 
    return size(eval.v1);
}


template <typename Scalar, unsigned long Unroll, typename Vector1, 
	  typename Vector2, typename ConjOpt, typename Assign>
dot_index_evaluator<Scalar, Vector1, Vector2, ConjOpt, Assign>
inline index_evaluator(lazy_assign<Scalar, mtl::dot_class<Unroll, Vector1, Vector2, ConjOpt>, Assign>& lazy)
{
    return dot_index_evaluator<Scalar, Vector1, Vector2, ConjOpt, Assign>(lazy.first, lazy.second.v1, lazy.second.v2);
}

template <typename VectorOut, typename Matrix, typename VectorIn, typename Assign>
struct row_mat_cvec_index_evaluator
{
    BOOST_STATIC_ASSERT((mtl::traits::is_row_major<Matrix>::value));
    typedef typename mtl::Collection<VectorOut>::value_type        value_type;
    typedef typename mtl::Collection<Matrix>::size_type            size_type; 

    row_mat_cvec_index_evaluator(VectorOut& w, const Matrix& A, const VectorIn& v) : w(w), A(A), v(v) {}

    template <unsigned Offset>
    void at(size_type i, mtl::tag::sparse)
    {
	value_type tmp(math::zero(w[i+Offset]));
	const size_type cj0= A.ref_major()[i+Offset], cj1= A.ref_major()[i+Offset+1];
	for (size_type j= cj0; j != cj1; ++j)
	    tmp+= A.data[j] * v[A.ref_minor()[j]];
	Assign::first_update(w[i+Offset], tmp);
    }

    template <unsigned Offset>
    void at(size_type i, mtl::tag::dense)
    {
	value_type tmp(math::zero(w[i+Offset]));
	for (size_type j= 0; j < num_cols(A); j++) 
	    tmp+= A[i][j] * v[j];
	Assign::first_update(w[i+Offset], tmp);
    }

    template <unsigned Offset>
    void at(size_type i)
    { 
	at<Offset>(i, typename mtl::traits::category<Matrix>::type());
    }

    void operator()(size_type i) { at<0>(i); }
    void operator[](size_type i) { at<0>(i); }

    VectorOut&      w;
    const Matrix&   A;
    const VectorIn& v;
};

template <typename VectorOut, typename Matrix, typename VectorIn, typename Assign>
inline std::size_t size(const row_mat_cvec_index_evaluator<VectorOut, Matrix, VectorIn, Assign>& eval)
{
    return size(eval.w);
}

template <typename VectorOut, typename Matrix, typename VectorIn, typename Assign>
row_mat_cvec_index_evaluator<VectorOut, Matrix, VectorIn, Assign>
inline index_evaluator(lazy_assign<VectorOut, mtl::mat_cvec_times_expr<Matrix, VectorIn>, Assign>& lazy)
{
    return row_mat_cvec_index_evaluator<VectorOut, Matrix, VectorIn, Assign>(lazy.first, lazy.second.first, lazy.second.second);
}

template <typename T, typename U> struct fused_expr;
template <typename T, typename U> void inline evaluate_lazy(fused_expr<T, U>& expr);


template <typename T, typename U>
struct fused_expr
{
    template <typename TT, typename UU, typename Assign>
    void check(lazy_assign<TT, UU, Assign>& )
    {
	bool vec_scal= boost::mpl::and_<mtl::traits::is_vector<TT>, mtl::traits::is_scalar<UU> >::value;
	bool vec_vec= boost::mpl::and_<mtl::traits::is_vector<TT>, mtl::traits::is_vector<UU> >::value;
	bool scal_red= boost::mpl::and_<mtl::traits::is_scalar<TT>, is_vector_reduction<UU> >::value;
	
	bool ia= boost::mpl::or_<
	    boost::mpl::and_<mtl::traits::is_vector<TT>, mtl::traits::is_scalar<UU> >,
	    boost::mpl::and_<mtl::traits::is_vector<TT>, mtl::traits::is_vector<UU> >,
	    boost::mpl::and_<mtl::traits::is_scalar<TT>, is_vector_reduction<UU> >
	    >::value;
    }


    fused_expr(T& first, U& second) : first(first), second(second) 
    {
	// check(first); check(second);
	//index_evaluatable<T> it= "";
	//index_evaluatable<U> iu= "";
    }
 
    ~fused_expr() { eval(index_evaluatable<T>(), index_evaluatable<U>()); }

    template <typename TT, typename UU>
    void eval_loop(TT first_eval, UU second_eval)
    {	
	MTL_DEBUG_THROW_IF(/*mtl::*/  size(first_eval) != /*mtl::*/  size(second_eval), mtl::incompatible_size());	

#ifdef MTL_LAZY_LOOP_WO_UNROLL
	for (std::size_t i= 0, s= size(first_eval); i < s; i++) {
	    first_eval(i); second_eval(i);
	}	
#else
	std::size_t s= size(first_eval), sb= s >> 2 << 2;

	for (std::size_t i= 0; i < sb; i+= 4) {
	    first_eval.template at<0>(i); second_eval.template at<0>(i);
	    first_eval.template at<1>(i); second_eval.template at<1>(i);
	    first_eval.template at<2>(i); second_eval.template at<2>(i);
	    first_eval.template at<3>(i); second_eval.template at<3>(i);
	}

	for (std::size_t i= sb; i < s; i++) {
	    first_eval(i); second_eval(i);
	}
#endif
    }

    void eval(boost::mpl::true_, boost::mpl::true_)
    {
	cout << "Now I really fuse!\n";
	eval_loop(index_evaluator(first), index_evaluator(second)); 
    }

    template <bool B1, bool B2>
    void eval(boost::mpl::bool_<B1>, boost::mpl::bool_<B2>)
    { evaluate_lazy(first); evaluate_lazy(second); }

    T& first;
    U& second;
};

template <typename T, typename U>
struct is_lazy<fused_expr<T, U> > 
  : boost::mpl::and_<is_lazy<T>, is_lazy<U> > 
{};

template <typename T, typename U>
struct index_evaluatable<fused_expr<T, U> > 
  : boost::mpl::and_<index_evaluatable<T>, index_evaluatable<U> > 
{};

template <typename T, typename U>
void inline evaluate_lazy(fused_expr<T, U>& expr) 
{ evaluate_lazy(expr.first); evaluate_lazy(expr.second); }


template <typename T, typename U>
struct fused_index_evaluator
{
    fused_index_evaluator(T& first, U& second) 
      : first(index_evaluator(first)), second(index_evaluator(second)) {}

    template <unsigned Offset>
    void at(std::size_t i) 
    { first.at<Offset>(i); second.at<Offset>(i); }

    void operator() (std::size_t i) { at<0>(i); }
    void operator[] (std::size_t i) { at<0>(i); }

    typename evaluator_type<T>::type first;
    typename evaluator_type<U>::type second;
};

template <typename T, typename U>
inline size_t size(const fused_index_evaluator<T, U>& expr) { return size(expr.first); }

template <typename T, typename U>
struct evaluator_type<fused_expr<T, U> >
{
    typedef fused_index_evaluator<T, U> type;
};

template <typename T, typename U>
inline fused_index_evaluator<T, U> index_evaluator(fused_expr<T, U>& expr)
{  return fused_index_evaluator<T, U>(expr.first, expr.second); }


template <typename T, typename U>
typename boost::enable_if<boost::mpl::and_<is_lazy<T>, is_lazy<U> >, fused_expr<T, U> >::type
operator||(const T& x, const U& y)
{
    return fused_expr<T, U>(const_cast<T&>(x), const_cast<U&>(y));
}

template <typename T, typename U>
typename boost::enable_if<boost::mpl::and_<is_lazy<T>, is_lazy<U> >, fused_expr<T, U> >::type
fuse(const T& x, const U& y)
{
    return fused_expr<T, U>(const_cast<T&>(x), const_cast<U&>(y));
}


int main(int, char**) 
{
    double                d, rho, alpha= 7.8, beta, gamma;
    const double          cd= 2.6;
    std::complex<double>  z;

    mtl::dense_vector<double> v(6, 1.0), w(6), r(6, 6.0), q(6, 2.0), x(6);
    mtl::dense2D<double>      A(6, 6);
    A= 2.0;
    mtl::compressed2D<double>      B(6, 6);
    B= 2.0;

    (lazy(w)= A * v) || (lazy(d) = lazy_dot(w, v));
    // fuse(lazy(w)= A * v, lazy(d) = lazy_dot(w, v));
    // d= with_reduction(lazy(w)= A * v, lazy_dot(w, v));
    cout << "w = " << w << ", d (12?)= " << d << "\n";

    (lazy(w)= B * v) || (lazy(d) = lazy_dot(w, v));
    // fuse(lazy(w)= A * v, lazy(d) = lazy_dot(w, v));
    // d= with_reduction(lazy(w)= A * v, lazy_dot(w, v));
    cout << "w = " << w << ", d (12?)= " << d << "\n";

    (lazy(r)-= alpha * q) || (lazy(rho)= lazy_unary_dot(r)); 
    //fuse( lazy(r)-= alpha * q, lazy(rho)= lazy_unary_dot(r) ); 
    // lazy(r)-= alpha * q, lazy(rho)= lazy_unary_dot(r);
    cout << "r = " << r << ", rho (552.96?) = " << rho << "\n";

    (lazy(x)= 7.0) || (lazy(beta)= lazy_unary_dot(x)); 
    cout << "x = " << x << ", beta (294?) = " << beta << "\n";
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_one_norm(x)); 
    cout << "x = " << x << ", beta (42?) = " << beta << "\n";
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_two_norm(x)); 
    cout << "x = " << x << ", beta (17.1464?) = " << beta << "\n";
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_infinity_norm(x)); 
    cout << "x = " << x << ", beta (7?) = " << beta << "\n";
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_sum(x)); 
    cout << "x = " << x << ", beta (42?) = " << beta << "\n";
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_product(x)); 
    cout << "x = " << x << ", beta (117649?) = " << beta << "\n";
    
    (lazy(x)= 2.0) || (lazy(gamma)= lazy_dot(r, x)); 
    cout << "x = " << x << ", gamma (-115.2?) = " << gamma << "\n";
    
    (lazy(r)= alpha * q) || (lazy(rho)= lazy_dot(r, q)); 
    cout << "r = " << r << ", rho (187.2?) = " << rho << "\n";

    (lazy(r)= alpha * q) || (lazy(v)= 8.6 * q) || (lazy(x)= 2.2 * q); 
    //fuse( lazy(r)-= alpha * q, lazy(rho)= lazy_unary_dot(r) ); 
    // lazy(r)-= alpha * q, lazy(rho)= lazy_unary_dot(r);
    cout << "r = " << r << ", v (17.2?) = " << v << "\n";


    


    return 0;
}
