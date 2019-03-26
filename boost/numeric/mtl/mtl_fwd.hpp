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

#ifndef MTL_MTL_FWD_INCLUDE
#define MTL_MTL_FWD_INCLUDE

/// Main name space for %Matrix Template Library
namespace mtl {
    
    template <typename T, typename U, typename Assign> struct lazy_assign;
    template <typename VectorOut, typename Matrix, typename VectorIn, typename Assign> struct row_mat_cvec_index_evaluator;

    /// Namespace for tags used for concept-free dispatching
    namespace tag {
	struct row_major;
	struct col_major;

	struct scalar;
	struct vector;
	struct matrix;

	/// Namespace for constant iterator tags
	namespace const_iter {}

	/// Namespace for iterator tags
	namespace iter {}
    }
    using tag::row_major;
    using tag::col_major;

    namespace index {
	struct c_index;
	struct f_index;
    }

    /// Namespace for compile-time parameters, e.g. %matrix dimensions
    namespace fixed {
	template <std::size_t Rows, std::size_t Cols> struct dimensions;
    }

    /// Namespace for run-time parameters, e.g. %matrix dimensions
    namespace non_fixed {
	struct dimensions;
    }

    /// Namespace for matrices and views and operations exclusively on matrices
    namespace mat {

	template <typename Orientation, typename Index, typename Dimensions, bool OnStack, typename SizeType> struct parameters;

        template <typename Value, typename Parameters> class dense2D;

        template <typename Value, typename Parameters> 
        typename dense2D<Value, Parameters>::size_type num_cols(const dense2D<Value, Parameters>& matrix);
        template <typename Value, typename Parameters> 
        typename dense2D<Value, Parameters>::size_type num_rows(const dense2D<Value, Parameters>& matrix);
        template <typename Value, typename Parameters> 
        typename dense2D<Value, Parameters>::size_type size(const dense2D<Value, Parameters>& matrix);

        template <typename Value, std::size_t Mask, typename Parameters> class morton_dense;

#if !defined(_MSC_VER) || _MSC_VER != 1400 // Trouble in MSVC 2005
        template <typename Value, std::size_t Mask, typename Parameters>
        typename morton_dense<Value, Mask, Parameters>::size_type num_cols(const morton_dense<Value, Mask, Parameters>& matrix);
        template <typename Value, std::size_t Mask, typename Parameters>
        typename morton_dense<Value, Mask, Parameters>::size_type num_rows(const morton_dense<Value, Mask, Parameters>& matrix);
        template <typename Value, std::size_t Mask, typename Parameters>
        typename morton_dense<Value, Mask, Parameters>::size_type size(const morton_dense<Value, Mask, Parameters>& matrix);
#endif

        template <typename Value, typename Parameters> class compressed2D;

        template <typename Value, typename Parameters> 
        typename compressed2D<Value, Parameters>::size_type num_cols(const compressed2D<Value, Parameters>& matrix);
        template <typename Value, typename Parameters> 
        typename compressed2D<Value, Parameters>::size_type num_rows(const compressed2D<Value, Parameters>& matrix);
        template <typename Value, typename Parameters> 
        // typename compressed2D<Value, Parameters>::size_type 
	std::size_t
	size(const compressed2D<Value, Parameters>& matrix);

        template <typename Value, typename Parameters, typename Updater> struct compressed2D_inserter;

	template <typename T, typename Parameters> class coordinate2D;
	template <typename Matrix, typename Updater> struct coordinate2D_inserter;

	template <typename T, typename Parameters> class sparse_banded;
	template <typename T, typename Parameters, typename Updater> struct sparse_banded_inserter;

	template <typename Value, typename Parameters> class ell_matrix;
	template <typename Value, typename Parameters, typename Updater> struct ell_matrix_inserter;

	template <typename Matrix, typename Updater> struct inserter;
	template <typename BaseInserter> class shifted_inserter;  

	template <typename Vector> class multi_vector;
	template <typename Vector> class multi_vector_range;
	template <typename Value> class element;
	template <typename Value> class element_structure;
	template <typename Functor> class implicit_dense;  
	template <typename Value> class ones_functor;
	template <typename Value> class ones_matrix;
	template <typename Value> class hilbert_functor;
	template <typename Value> class hilbert_matrix;
	template <typename Vector1, typename Vector2> class outer_product_functor;
	template <typename Vector1, typename Vector2> class outer_product_matrix;

        struct identity2D;
	std::size_t size(const identity2D& A);
	std::size_t num_rows(const identity2D& A);
	std::size_t num_cols(const identity2D& A);

        template <typename Matrix> struct transposed_view;
        
	template <typename Matrix> struct mat_expr;
	template <typename Matrix> struct dmat_expr;
	template <typename Matrix> struct smat_expr;
	template <typename M1, typename M2, typename SFunctor> struct mat_mat_op_expr;
	template <typename M1, typename M2> struct mat_mat_plus_expr;
	template <typename M1, typename M2> struct mv_mv_plus_expr;
	template <typename M1, typename M2> struct mat_mat_minus_expr;
	template <typename M1, typename M2> struct mv_mv_minus_expr;
	template <typename M1, typename M2> struct mat_mat_ele_times_expr;
	template <typename M1, typename M2> struct mat_mat_times_expr;
	template <typename M1, typename M2> struct mat_mat_asgn_expr;

	template <typename Matrix> struct mat_expr;
	template <typename Functor, typename Matrix> struct map_view;
	template <typename Scaling, typename Matrix> struct scaled_view;
	template <typename Matrix, typename RScaling> struct rscaled_view; // added by Hui Li
	template <typename Matrix, typename Divisor> struct divide_by_view; // added by Hui Li
	template <typename Matrix>  struct conj_view;
	template <typename Matrix>  struct negate_view;
	template <typename Matrix>  struct imag_view;
	template <typename Matrix>  struct real_view;
	template <typename Matrix>  struct hermitian_view;
	template <typename Matrix>  struct banded_view;
	template <typename Matrix> struct exp_view;
	template <typename Matrix> struct indirect;

	template <typename Matrix> class lu_solver;

	template <typename Matrix> std::size_t size(const banded_view<Matrix>&);
	template <class Matrix> std::size_t size(const transposed_view<Matrix>&);
	template <typename Functor, typename Matrix> std::size_t size(const map_view<Functor, Matrix>&);
	template <class Matrix> std::size_t size(const hermitian_view<Matrix>& );
    }

    //using mat::dense2D;
    //using mat::morton_dense;
    //using mat::compressed2D;
    //using mat::coordinate2D;
    //using mat::multi_vector;
    //using mat::transposed_view;
    
    template <typename E1, typename E2> struct mat_cvec_times_expr;

    /// Namespace for vectors and views and %operations exclusively on vectors
    namespace vec {
	template <typename Vector> struct vec_expr;
	template <typename Value, typename Parameters> class dense_vector;
	template <typename Value, typename Parameters> class strided_vector_ref;
	template <typename Value, typename Parameters> class sparse_vector;
	template <typename Functor, typename Vector> struct map_view;
	template <typename Vector, typename Exponent> struct pow_by_view;
	template <typename Vector>  struct conj_view;
	template <typename Vector>  struct real_view;
	template <typename Vector>  struct imag_view;
	template <typename Vector>  struct negate_view;
	template <typename Vector>  struct abs_view;
	template <typename Vector>  struct acos_view;
	template <typename Vector>  struct acosh_view;
	template <typename Vector>  struct asin_view;
	template <typename Vector>  struct asinh_view;
	template <typename Vector>  struct atan_view;
	template <typename Vector>  struct atanh_view;
	template <typename Vector>  struct cos_view;
	template <typename Vector>  struct cosh_view;
	template <typename Vector>  struct sin_view;
	template <typename Vector>  struct sinh_view;
	template <typename Vector>  struct tan_view;
	template <typename Vector>  struct tanh_view;
	template <typename Vector>  struct ceil_view;
	template <typename Vector>  struct floor_view;
	template <typename Vector>  struct log_view;
	template <typename Vector>  struct log10_view;
	template <typename Vector>  struct exp_view;
	template <typename Vector>  struct exp10_view;
	template <typename Vector>  struct sqrt_view;
	template <typename Vector>  struct rsqrt_view;
	template <typename Vector>  struct signum_view;
  
# ifdef MTL_WITH_MATH_ELEVEN    
	template <typename Vector>  struct round_view;
	template <typename Vector>  struct trunc_view;
	template <typename Vector>  struct log2_view;
	template <typename Vector>  struct exp2_view;
	template <typename Vector>  struct erf_view;
	template <typename Vector>  struct erfc_view;  
# endif
  
	template <typename Scaling, typename Vector> struct scaled_view;
	template <typename Vector, typename RScaling> struct rscaled_view; // added by Hui Li
	template <typename Vector, typename Divisor> struct divide_by_view; // added by Hui Li
	template <class E1, class E2, typename SFunctor> struct vec_vec_op_expr;
	template <class E1, class E2, typename SFunctor> struct vec_vec_pmop_expr;
	template <class E1, class E2, typename SFunctor> struct vec_vec_aop_expr;
	template <class E1, class E2, typename SFunctor> struct vec_vec_ele_prod_expr;
	template <class E1, class E2, typename SFunctor> struct vec_scal_aop_expr;
	template <class E1, class E2> struct vec_vec_asgn_expr;
	template <class E1, class E2> struct vec_vec_plus_asgn_expr;
	template <class E1, class E2> struct vec_vec_minus_asgn_expr;
	template <class E1, class E2> struct vec_vec_times_asgn_expr; // is this really used???
	//  template <class E1, class E2> struct vec_vec_div_asgn_expr; // is this really used???
	template <class E1, class E2> struct vec_scal_times_asgn_expr;
	template <class E1, class E2> struct vec_scal_div_asgn_expr; // added by Hui Li
	template <class E1, class E2> struct vec_scal_asgn_expr;
	template <typename E1, typename E2> struct rvec_mat_times_expr;

	template <typename Vector> struct vec_const_ref_expr;
	template <unsigned BSize, typename Vector> class unrolled1;  

	template <typename Scalar, typename Vector, typename Functor, typename Assign> struct reduction_index_evaluator;
	template <typename Scalar, typename Vector1, typename Vector2, typename ConjOpt, typename Assign> struct dot_index_evaluator;
	template <unsigned long Unroll, typename Vector1, typename Vector2, typename ConjOpt> struct dot_class;

	template <typename Vector, typename Functor> struct lazy_reduction;

	template <typename Matrix, typename VectorIn> struct mat_cvec_multiplier;

	template <typename Value, typename Parameters, typename Value2>
	inline void fill(dense_vector<Value, Parameters>& vector, const Value2& value);
  
	template <typename Value, typename Parameters>
	typename dense_vector<Value, Parameters>::size_type
	inline size(const dense_vector<Value, Parameters>& vector);

	template <typename Value, typename Parameters>
	typename dense_vector<Value, Parameters>::size_type
	inline num_rows(const dense_vector<Value, Parameters>& vector);
  
	template <typename Value, typename Parameters>
	typename dense_vector<Value, Parameters>::size_type
	inline num_cols(const dense_vector<Value, Parameters>& vector);

	template <typename Value, typename Parameters>
	std::size_t size(const strided_vector_ref<Value, Parameters>& v);

	template <typename Functor, typename Vector> 
	std::size_t size(const map_view<Functor, Vector>& v);

	template <typename E1, typename E2, typename SFunctor>
	std::size_t size(const vec_vec_op_expr<E1, E2, SFunctor>& v);

	template <typename E1, typename E2, typename SFunctor>
	std::size_t size(const vec_vec_aop_expr<E1, E2, SFunctor>& v);

	template <typename E1, typename E2, typename SFunctor>
	std::size_t size(const vec_vec_pmop_expr<E1, E2, SFunctor>& v);

	template <typename E1, typename E2, typename SFunctor>
	inline std::size_t size(const vec_scal_aop_expr<E1, E2, SFunctor>& v);

	template <unsigned BSize, typename Vector>
	inline std::size_t size(const unrolled1<BSize, Vector>& v);

	template <typename E1, typename E2>
	std::size_t inline size(const rvec_mat_times_expr<E1, E2>& x);

	template <typename E1, typename E2>
	std::size_t inline size(const mat_cvec_times_expr<E1, E2>& x);

	template <typename Matrix, typename VectorIn>
	std::size_t size(const mat_cvec_multiplier<Matrix, VectorIn>& m);

	/// Namespace for fixed vector dimension types
	namespace fixed {
	    template <std::size_t Size> struct dimension;
	}
	/// Namespace for non-fixed vector dimension types, i.e. size dynamically determined at run-time
	namespace non_fixed {
	    struct dimension;
	}
      
    }

    //using dense_vector;

    // Export free vector functions into mtl namespace
    // It is also needed to handle STL vectors in MTL
    using vec::fill;
    using vec::size;
    using vec::num_rows;
    using vec::num_cols;

    /// Namespace for %operations (if not defined in mtl)
    namespace operations {
	template <typename T> struct update_store;
    }


    namespace vec {

	template <typename Vector, typename Updater = mtl::operations::update_store<typename Vector::value_type> > struct inserter;
	template <typename Vector, typename Size> struct update_proxy;
    }





    /// Namespace for type %traits
    namespace traits {
	template <typename Value> struct category;
	template <typename Value> struct algebraic_category;

	template <typename Collection> struct value;
	template <typename Collection> struct const_value;
	template <typename Collection> struct row;
	template <typename Collection> struct col;
	template <class Matrix> struct offset;

	template <class Vector> struct index;
	template <typename Tag, typename Collection>  struct range_generator;

	template <typename T> struct eval_dense;

        template <typename Matrix> struct transposed_orientation;

	template <typename Collection> struct root;

	// for internal implementations
	namespace detail {
	    // needed collection.hpp (at least)
	    template <typename Collection, typename Cursor, typename Complexity> struct dense_element_range_generator;
	    template <typename Matrix, typename Cursor, typename Complexity>     struct all_offsets_range_generator;
	    template <typename Matrix, typename Tag, int Level = 2>              struct sub_matrix_cursor;
	    template <typename Matrix>                                           struct matrix_element_key;
	    template <typename Matrix, int pos>                                  struct matrix_element_cursor;
	    template <typename Matrix, typename Complexity, int Level = 2>       struct all_rows_range_generator;
	    template <typename Cursor>                                           struct all_cols_in_row_range_generator;
	    template <typename Matrix, typename Complexity, int Level = 2>       struct all_cols_range_generator;
	    template <typename Cursor>                                           struct all_rows_in_col_range_generator;
	    template <typename Collection, typename RangeGenerator>              struct referred_range_generator;
	}
    }

    template <typename Matrix> struct ColumnInMatrix;
    template <typename Matrix> struct RowInMatrix;

    template <class Tag, class Collection> typename traits::range_generator<Tag, Collection>::type 
    begin(Collection const& c);
    
    template <class Tag, class Collection> typename traits::range_generator<Tag, Collection>::type 
    end(Collection const& c);


    /// Namespace for functors with application operator and fully typed parameters
    namespace tfunctor {
	/// Functor for scaling matrices, vectors and ordinary scalars
	template <typename V1, typename V2, typename AlgebraicCategory = tag::scalar> struct scale;

	/// Functor for scaling matrices, vectors and ordinary scalars
	template <typename V1, typename V2, typename AlgebraicCategory = tag::scalar> struct rscale;

	/// Functor for scaling matrices, vectors and ordinary scalars
	template <typename V1, typename V2, typename AlgebraicCategory = tag::scalar> struct divide_by;
        
        template <typename Value1, typename Value2> struct pow_by;
    }

    /// Namespace for functors with static function apply and fully typed parameters
    namespace sfunctor {
	template <typename Value, typename AlgebraicCategory = tag::scalar> struct conj_aux;
	template <typename Value>                   struct conj;
	template <typename Value>                   struct imag;
	template <typename Value>                   struct real;
	template <typename Value>                   struct negate;
	template <typename Value1, typename Value2> struct plus;
	template <typename Value1, typename Value2> struct minus;
	template <typename Value1, typename Value2> struct times;
	template <typename Value1, typename Value2> struct divide;
	template <typename Value1, typename Value2> struct assign;
	template <typename Value1, typename Value2> struct plus_assign;
	template <typename Value1, typename Value2> struct minus_assign;
	template <typename Value1, typename Value2> struct times_assign;
	template <typename Value1, typename Value2> struct divide_assign;
	template <typename Value>                   struct identity;
	template <typename Value>                   struct abs;
	template <typename Value>                   struct exp;
	template <typename Value>                   struct sqrt;
	template <typename Value>                   struct square;
	template <typename F, typename G>           struct compose;
	template <typename F, typename G>           struct compose_first;
	template <typename F, typename G>           struct compose_second;
	template <typename F, typename G, typename H> struct compose_both;
	template <typename F, typename G>           struct compose_binary;
    }

    // Namespace documentations

    /// Namespace for static assignment functors
    namespace assign {}

    /// Namespace for complexity classes
    namespace complexity_classes {}

    /// Namespace for %operations (if not defined in mtl)
    namespace operations {}

    /// Namespace for recursive operations and types with recursive memory layout
    namespace recursion {}

    /// Namespace for %utilities
    namespace utility {}

    /// Namespace for implementations using recursators
    namespace wrec {}

    namespace mat {
	template <typename Matrix, typename ValueType, typename SizeType> struct crtp_matrix_assign;
	template <typename Matrix, typename ValueType, typename SizeType> struct const_crtp_matrix_bracket;
	template <typename Matrix, typename ValueType, typename SizeType> struct crtp_matrix_bracket;
	template <typename Matrix, typename ValueType, typename SizeType> struct const_crtp_base_matrix;
	template <typename Matrix, typename ValueType, typename SizeType> struct crtp_base_matrix;
	template <typename Matrix> struct const_crtp_matrix_range_bracket;
    }

    namespace detail {
	template <typename Value, bool OnStack, unsigned Size= 0> struct contiguous_memory_block;
	template <typename Matrix, typename Updater> struct trivial_inserter;
	template <typename Collection> struct with_format_t;
    }

    /// Free function defined for all matrix and vector types
    template <typename Collection> void swap(Collection& c1, Collection& c2);

    /// User registration that class has a clone constructor, otherwise use regular copy constructor.
    template<typename T> struct is_clonable;

    /// Helper type to define constructors that always copy
    struct clone_ctor {};

    /// Helper type to define constructors that refer to another object instead of copying it (e.g. transposed vectors)
    struct refer_ctor {};

    class irange;
    class iset;
    class srange;

    template <typename T, typename U> struct fused_expr;
    namespace vec {
	template <typename T, typename U> struct fused_index_evaluator;
	// template <typename T, typename U> size_t size(const fused_index_evaluator<T, U>&); // not needed currently
    }

    /// Namespace for I/O operations
    namespace io {
	class matrix_market_istream;
	class matrix_market_ostream;

	template <typename MatrixIStream, typename MatrixOStream> class matrix_file;
	typedef matrix_file<matrix_market_istream, matrix_market_ostream> matrix_market;
    }

    // Multiplication functors
    template <typename Assign, typename Backup> struct gen_cursor_dmat_dmat_mult_t;
    template <typename Assign, typename Backup> struct gen_dmat_dmat_mult_t;
    template <unsigned long Tiling1, unsigned long Tiling2, typename Assign, typename Backup> struct gen_tiling_dmat_dmat_mult_t;
    template <typename Assign, typename Backup> struct gen_tiling_44_dmat_dmat_mult_t;
    template <typename Assign, typename Backup> struct gen_tiling_22_dmat_dmat_mult_t;
    template <typename BaseMult, typename BaseTest, typename Assign, typename Backup> struct gen_recursive_dmat_dmat_mult_t;
    template <typename Assign, typename Backup> struct gen_platform_dmat_dmat_mult_t;
    template <typename Assign, typename Backup> struct gen_blas_dmat_dmat_mult_t;
    template <std::size_t SizeLimit, typename FunctorSmall, typename FunctorLarge> struct size_switch_dmat_dmat_mult_t;
    template <bool IsStatic, typename FunctorStatic, typename FunctorDynamic> struct static_switch_dmat_dmat_mult_t;
    template <typename Assign, typename Backup> struct fully_unroll_fixed_size_dmat_dmat_mult_t;

} // namespace mtl

#endif // MTL_MTL_FWD_INCLUDE


