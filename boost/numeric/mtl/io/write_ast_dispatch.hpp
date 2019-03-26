// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_IO_WRITE_AST_DISPATCH_INCLUDE
#define MTL_IO_WRITE_AST_DISPATCH_INCLUDE

#include <string>
#include <sstream>
#include <fstream>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/io/functor_symbol.hpp>
#include <boost/numeric/mtl/vector/vec_vec_aop_expr.hpp>

namespace mtl { namespace io {

template <typename Value, typename Parameters> 
void write_ast_dispatch(const mtl::dense_vector<Value, Parameters>& v, std::string s, std::ofstream& f);
template <typename E1, typename E2, typename SFunctor>
void write_ast_dispatch(const mtl::vec_vec_aop_expr<E1, E2, SFunctor>& expr, std::string s, std::ofstream& f);
template <class E1, class E2, typename SFunctor>
void write_ast_dispatch(const mtl::vec_vec_pmop_expr<E1, E2, SFunctor>& expr, std::string s, std::ofstream& f);
template <typename Functor, typename Vector> 
void write_ast_dispatch(const mtl::map_view<Functor, Vector>& expr, std::string s, std::ofstream& f);
template <typename Scaling, typename Vector>
void write_ast_dispatch(const mtl::scaled_view<Scaling, Vector>& expr, std::string s, std::ofstream& f);
template <typename Matrix, typename Vector> 
void write_ast_dispatch(const mtl::mat_cvec_times_expr<Matrix, Vector>& expr, std::string s, std::ofstream& f);
template <typename Expr>
void write_ast_dispatch(const mtl::operation::compute_summand<Expr>& expr, std::string s, std::ofstream& f);
template <typename Matrix, typename Vector> 
void write_ast_dispatch(const mtl::operation::compute_summand<mtl::mat_cvec_times_expr<Matrix, Vector> >& expr, std::string s, std::ofstream& f);
template <typename Vector1, typename Vector2> 
void write_ast_dispatch(const mtl::mat::outer_product_matrix<Vector1, Vector2>& expr, std::string s, std::ofstream& f);



template <typename Value>
typename boost::enable_if_c<boost::is_floating_point<Value>::value || boost::is_integral<Value>::value>::type
write_ast_dispatch(const Value& v, std::string s, std::ofstream& f)
{   
    f << "  " << s << "[shape=box,label=\"scalar\\n" << v << "\"]\n";
}


template <typename Value, typename Parameters> 
void write_ast_dispatch(const mtl::dense_vector<Value, Parameters>& v, std::string s, std::ofstream& f)
{   
    f << "  " << s << "[shape=box,label=\"vector\\n" << &v << "\"]\n";
}

template <typename E1, typename E2, typename SFunctor>
void write_ast_dispatch(const mtl::vec_vec_aop_expr<E1, E2, SFunctor>& expr, std::string s, std::ofstream& f)
{
    f << "  " << s << "[label=\"" << functor_symbol(SFunctor()) << "\"]\n"; 
    std::string target= s + "t", source= s + "s";
    write_ast_dispatch(expr.first_argument(), target, f);
    write_ast_dispatch(expr.second_argument(), source, f);
    f << "  " << s << "->" << target << '\n';
    f << "  " << s << "->" << source << '\n';
}

template <class E1, class E2, typename SFunctor>
void write_ast_dispatch(const mtl::vec_vec_pmop_expr<E1, E2, SFunctor>& expr, std::string s, std::ofstream& f)
{
    f << "  " << s << "[label=\"" << functor_symbol(SFunctor()) << "\"]\n"; 
    std::string first= s + "f", second= s + "s";
    write_ast_dispatch(expr.first_argument(), first, f);
    write_ast_dispatch(expr.second_argument(), second, f);
    f << "  " << s << "->" << first << '\n';
    f << "  " << s << "->" << second << '\n';
}

template <typename Scaling, typename Vector>
void write_ast_dispatch(const mtl::scaled_view<Scaling, Vector>& expr, std::string s, std::ofstream& f)
{
    f << "  " << s << "[label=\"scaled_view\"]\n"; 
    std::string functor= s + "f", ref= s + "r";
    
    write_ast_dispatch(expr.functor.value, functor, f);
    write_ast_dispatch(expr.ref, ref, f);

    f << "  " << s << "->" << functor << '\n';
    f << "  " << s << "->" << ref << '\n';
}

template <typename Value1, typename Value2>
void write_ast_dispatch(const mtl::tfunctor::scale<Value1, Value2, mtl::tag::scalar>& expr, std::string s, std::ofstream& f)
{
    f << "  " << s << "[label=\"scale\"]\n"; 
    std::string factor= s + "f", ref= s + "r";
    write_ast_dispatch(expr.value(), factor, f);
    f << "  " << ref << "[label=\".\"]\n"; 

    f << "  " << s << "->" << factor << '\n';
    f << "  " << s << "->" << ref << '\n';
}

template <typename Functor, typename Vector> 
void write_ast_dispatch(const mtl::map_view<Functor, Vector>& expr, std::string s, std::ofstream& f)
{
    f << "  " << s << "[label=\"map\"]\n"; 
    std::string functor= s + "f", ref= s + "r";
    
    write_ast_dispatch(expr.functor, functor, f);
    write_ast_dispatch(expr.ref, ref, f);

    f << "  " << s << "->" << functor << '\n';
    f << "  " << s << "->" << ref << '\n';
}

// template <typename Matrix, typename Vector> 
// void write_ast_dispatch(const mtl::mat_cvec_times_expr<Matrix, Vector>& expr, std::string s, std::ofstream& f)
// {
// }

template <typename Expr>
void write_ast_dispatch(const mtl::operation::compute_summand<Expr>& expr, std::string s, std::ofstream& f)
{
    write_ast_dispatch(expr.value, s, f);
}

template <typename Matrix, typename Vector> 
void write_ast_dispatch(const mtl::operation::compute_summand<mtl::mat_cvec_times_expr<Matrix, Vector> >& expr, std::string s, std::ofstream& f)
{
# ifdef NDEBUG
    write_ast_dispatch(expr.value, s, f);
# else
    f << "  " << s << "[label=\"*\"]\n"; 
    std::string first= s + "f", second= s + "s";
    write_ast_dispatch(expr.first, first, f);
    write_ast_dispatch(expr.second, second, f);
    f << "  " << s << "->" << first << '\n';
    f << "  " << s << "->" << second << '\n';
# endif 
}

template <typename Vector1, typename Vector2> 
void write_ast_dispatch(const mtl::mat::outer_product_matrix<Vector1, Vector2>& expr, std::string s, std::ofstream& f)
{
    f << "  " << s << "[label=\"outer product\"]\n"; 
    std::string first= s + "f", second= s + "s";
    write_ast_dispatch(expr.v1(), first, f);
    write_ast_dispatch(expr.v2(), second, f);
    f << "  " << s << "->" << first << '\n';
    f << "  " << s << "->" << second << '\n';
}

}} // namespace mtl::io

#endif // MTL_IO_WRITE_AST_DISPATCH_INCLUDE
