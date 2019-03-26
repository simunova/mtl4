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

#ifndef MTL_MATRIX_UMFPACK_SOLVE_INCLUDE
#define MTL_MATRIX_UMFPACK_SOLVE_INCLUDE

#ifdef MTL_HAS_UMFPACK

#include <iostream>


#include <cassert>
#include <algorithm>
#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/make_copy_or_reference.hpp>
#include <boost/numeric/mtl/operation/base_solver.hpp>
#include <boost/numeric/mtl/operation/merge_complex_vector.hpp>
#include <boost/numeric/mtl/operation/split_complex_vector.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

extern "C" {
#  include <umfpack.h>
}

namespace mtl { namespace mat {

    /// Namespace for Umfpack solver
    namespace umfpack {

	// conversion for value_type needed if not double or complex<double> (where possible)
	template <typename Value> struct value         	    {};
	template<> struct value<long double>           	    { typedef double               type; };
	template<> struct value<double>                	    { typedef double               type; };
	template<> struct value<float>                 	    { typedef double               type; };
	template<> struct value<std::complex<long double> > { typedef std::complex<double> type; };
	template<> struct value<std::complex<double> > 	    { typedef std::complex<double> type; };
	template<> struct value<std::complex<float> >  	    { typedef std::complex<double> type; };

	template <typename Value> struct use_long { static const bool value= sizeof(Value) > sizeof(int); };

	template <bool Larger> struct index_aux   { typedef int     type; };
#     ifdef UF_long
        template<> struct index_aux<true>         { typedef UF_long type; };
#     elif defined(SuiteSparse_long)
        template<> struct index_aux<true>         { typedef SuiteSparse_long type; };
#     else
#pragma message "cannot deduce the long umfpack type"
	template<> struct index_aux<true>         { typedef long type; };
#     endif

	template <typename Value> struct index 
          : index_aux<use_long<Value>::value> {};

	template <typename Matrix, typename Value, typename Orientation> 
	struct matrix_copy {};

	// If arbitrary compressed matrix -> copy
	template <typename Value, typename Parameters, typename Orientation>
	struct matrix_copy<compressed2D<Value, Parameters>, Value, Orientation>
	{
	    typedef typename value<Value>::type                      value_type;
	    typedef compressed2D<value_type, parameters<col_major> > matrix_type;
	    typedef compressed2D<Value, Parameters>                  in_matrix_type;

	    matrix_copy(const in_matrix_type& A) : matrix(A) {}
	    matrix_type matrix;
	};

	struct error : public domain_error
	{
	    error(const char *s, int code) : domain_error(s), code(code) {}
	    int code;
	};

	inline void check(int res, 
			  const char*
#ifndef MTL_ASSERT_FOR_THROW 
			  s
#endif
			  )
	{
#ifndef MTL_ASSERT_FOR_THROW 
	    MTL_THROW_IF(res != UMFPACK_OK, error(s, res));
#else
	    MTL_THROW_IF(res != UMFPACK_OK, error("", res));
#endif
	}

	/// Class for repeated Umfpack solutions
	/** Keeps symbolic and numeric preprocessing. Numeric part can be updated. 
	    Only defined for compressed2D<double> and compressed2D<complex<double> >. **/
	template <typename T> 
	class solver {
	  public:
	    /// Constructor referring to matrix \p A (not changed) and optionally Umfpack's strategy and alloc_init (look for the specializations)
	    // \ref solver<compressed2D<double, Parameters> > and \ref solver<compressed2D<std::complex<double>, Parameters> >)
	    explicit solver(const T& A) {}

	    /// Update numeric part, for matrices that kept the sparsity and changed the values
	    void update_numeric() {}

	    /// Update symbolic and numeric part
	    void update() {}

	    /// Solve system A*x == b with matrix passed in constructor
	    /** Please note that the order of b and x is different than in solve() !!! **/
	    template <typename VectorX, typename VectorB>
	    int operator()(VectorX& x, const VectorB& b) const {return 0;}

	    /// Solve system A*x == b with matrix passed in constructor
	    /** Please note that the order of b and x is different than in operator() !!! **/
	    template <typename VectorB, typename VectorX>
	    int solve(const VectorB& b, VectorX& x) const {return 0;}
	};

	/// Speciatization of solver for \ref mat::compressed2D with double values
	template <typename Parameters>
	class solver<compressed2D<double, Parameters> >
	  : public base_solver
	{
	    typedef double                                    value_type;
	    typedef compressed2D<value_type, Parameters>      matrix_type;
	    typedef typename matrix_type::size_type           size_type;
	    typedef typename index<size_type>::type           index_type;

	    static const bool copy_indices= sizeof(index_type) != sizeof(size_type),
		              long_indices= use_long<size_type>::value;
	    typedef boost::mpl::bool_<long_indices>           blong;
	    typedef boost::mpl::true_                         true_;
	    typedef boost::mpl::false_                        false_;

	    // typedef parameters<col_major>     Parameters;

	    void assign_pointers()
	    {
		if (copy_indices) {
		    if (Apc == 0) Apc= new index_type[n + 1]; 
		    if (my_nnz != A.nnz() && Aic) { delete[] Aic; Aic= 0; }
		    if (Aic == 0) Aic= new index_type[A.nnz()];
		    std::copy(A.address_major(), A.address_major() + n + 1, Apc);
		    std::copy(A.address_minor(), A.address_minor() + A.nnz(), Aic);
		    Ap= Apc;
		    Ai= Aic;
		} else {
		    Ap= reinterpret_cast<const index_type*>(A.address_major());
		    Ai= reinterpret_cast<const index_type*>(A.address_minor());
		}
		Ax= A.address_data();
	    }

	    void init_aux(true_)
	    {
	      check(umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &Symbolic, Control, Info), "Error in dl_symbolic");
	      check(umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info), "Error in dl_numeric");
	    }
	    
	    void init_aux(false_)
	    {
		check(umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, Control, Info), "Error in di_symbolic");
#if 0
		std::cout << "=== INFO of umfpack_*_symbolic ===\n";
		std::cout << "    UMFPACK_STATUS: " << (Info[UMFPACK_STATUS] == UMFPACK_OK ? "OK" : "ERROR") << "\n";
		std::cout << "    UMFPACK_NROW: " << Info[UMFPACK_NROW] << "\n";
		std::cout << "    UMFPACK_NCOL: " << Info[UMFPACK_NCOL] << "\n";
		std::cout << "    UMFPACK_NZ: " << Info[UMFPACK_NZ] << "\n";
		std::cout << "    UMFPACK_SIZE_OF_UNIT: " << Info[UMFPACK_SIZE_OF_UNIT] << "\n";
		std::cout << "    UMFPACK_NDENSE_ROW: " << Info[UMFPACK_NDENSE_ROW] << "\n";
		std::cout << "    UMFPACK_NEMPTY_ROW: " << Info[UMFPACK_NEMPTY_ROW] << "\n";
		std::cout << "    UMFPACK_NDENSE_COL: " << Info[UMFPACK_NDENSE_COL] << "\n";
		std::cout << "    UMFPACK_NEMPTY_COL: " << Info[UMFPACK_NEMPTY_COL] << "\n";
		std::cout << "    UMFPACK_SYMBOLIC_DEFRAG: " << Info[UMFPACK_SYMBOLIC_DEFRAG] << "\n";
		std::cout << "    UMFPACK_SYMBOLIC_PEAK_MEMORY: " << Info[UMFPACK_SYMBOLIC_PEAK_MEMORY] << "\n";
		std::cout << "    UMFPACK_SYMBOLIC_SIZE: " << Info[UMFPACK_SYMBOLIC_SIZE] << "\n";
		std::cout << "    UMFPACK_VARIABLE_PEAK_ESTIMATE: " << Info[UMFPACK_VARIABLE_PEAK_ESTIMATE]  << "\n";
		std::cout << "    UMFPACK_NUMERIC_SIZE_ESTIMATE: " << Info[UMFPACK_NUMERIC_SIZE_ESTIMATE] << "\n";
		std::cout << "    UMFPACK_PEAK_MEMORY_ESTIMATE: " << Info[UMFPACK_PEAK_MEMORY_ESTIMATE] << "\n";
		std::cout << "    UMFPACK_FLOPS_ESTIMATE: " << Info[UMFPACK_FLOPS_ESTIMATE] << "\n";
		std::cout << "    UMFPACK_LNZ_ESTIMATE: " << Info[UMFPACK_LNZ_ESTIMATE] << "\n";
		std::cout << "    UMFPACK_UNZ_ESTIMATE: " << Info[UMFPACK_UNZ_ESTIMATE] << "\n";
		std::cout << "    UMFPACK_MAX_FRONT_SIZE_ESTIMATE: " << Info[UMFPACK_MAX_FRONT_SIZE_ESTIMATE] << "\n";
		std::cout << "    UMFPACK_SYMBOLIC_TIME: " << Info[UMFPACK_SYMBOLIC_TIME] << "\n";
		std::cout << "    UMFPACK_SYMBOLIC_WALLTIME: " << Info[UMFPACK_SYMBOLIC_WALLTIME] << "\n";
		
		if (Info[UMFPACK_STRATEGY_USED] == UMFPACK_STRATEGY_SYMMETRIC)
		  std::cout << "    UMFPACK_STRATEGY_USED: SYMMETRIC\n";
		else {
		  if(Info[UMFPACK_STRATEGY_USED] == UMFPACK_STRATEGY_UNSYMMETRIC)
		    std::cout << "    UMFPACK_STRATEGY_USED: UNSYMMETRIC\n";
		  else {
		    if (Info[UMFPACK_STRATEGY_USED] == UMFPACK_STRATEGY_2BY2)
		      std::cout << "    UMFPACK_STRATEGY_USED: 2BY2\n";
		    else
		      std::cout << "    UMFPACK_STRATEGY_USED: UNKOWN STRATEGY " << Info[UMFPACK_STRATEGY_USED] << "\n";
		  }
		}
		  
		std::cout << "    UMFPACK_ORDERING_USED: " << Info[UMFPACK_ORDERING_USED] << "\n";
		std::cout << "    UMFPACK_QFIXED: " << Info[UMFPACK_QFIXED] << "\n";
		std::cout << "    UMFPACK_DIAG_PREFERRED: " << Info[UMFPACK_DIAG_PREFERRED] << "\n";
		std::cout << "    UMFPACK_ROW_SINGLETONS: " << Info[UMFPACK_ROW_SINGLETONS] << "\n";
		std::cout << "    UMFPACK_COL_SINGLETONS: " << Info[UMFPACK_COL_SINGLETONS] << "\n";
		std::cout << "    UMFPACK_PATTERN_SYMMETRY: " << Info[UMFPACK_PATTERN_SYMMETRY] << "\n";
		std::cout << "    UMFPACK_NZ_A_PLUS_AT: " << Info[UMFPACK_NZ_A_PLUS_AT] << "\n";
		std::cout << "    UMFPACK_NZDIAG: " << Info[UMFPACK_NZDIAG] << "\n";
		std::cout << "    UMFPACK_N2: " << Info[UMFPACK_N2] << "\n";
		std::cout << "    UMFPACK_S_SYMMETRIC: " << Info[UMFPACK_S_SYMMETRIC] << "\n";
		std::cout << "    UMFPACK_MAX_FRONT_NROWS_ESTIMATE: " << Info[UMFPACK_MAX_FRONT_NROWS_ESTIMATE] << "\n";
		std::cout << "    UMFPACK_MAX_FRONT_NCOLS_ESTIMATE: " << Info[UMFPACK_MAX_FRONT_NCOLS_ESTIMATE] << "\n";
		std::cout << "    UMFPACK_SYMMETRIC_LUNZ: " << Info[UMFPACK_SYMMETRIC_LUNZ] << "\n";
		std::cout << "    UMFPACK_SYMMETRIC_FLOPS: " << Info[UMFPACK_SYMMETRIC_FLOPS] << "\n";
		std::cout << "    UMFPACK_SYMMETRIC_NDENSE: " << Info[UMFPACK_SYMMETRIC_NDENSE] << "\n";
		std::cout << "    UMFPACK_SYMMETRIC_DMAX: " << Info[UMFPACK_SYMMETRIC_DMAX] << "\n";
#endif

		check(umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info), "Error in di_numeric");

#if 0
		std::cout << "=== INFO of umfpack_*_numeric ===\n";
		std::cout << "    UMFPACK_STATUS: " << (Info[UMFPACK_STATUS] == UMFPACK_OK ? "OK" : "ERROR") << "\n";
		std::cout << "    UMFPACK_VARIABLE_PEAK: " << Info[UMFPACK_VARIABLE_PEAK]  << "\n";
		std::cout << "    UMFPACK_PEAK_MEMORY: " << Info[UMFPACK_PEAK_MEMORY] << "\n";
		std::cout << "    UMFPACK_FLOPS: " << Info[UMFPACK_FLOPS] << "\n";
		std::cout << "    UMFPACK_LNZ: " << Info[UMFPACK_LNZ] << "\n";
		std::cout << "    UMFPACK_UNZ: " << Info[UMFPACK_UNZ] << "\n";
		std::cout << "    UMFPACK_NUMERIC_DEFRAG: " << Info[UMFPACK_NUMERIC_DEFRAG] << "\n";
		std::cout << "    UMFPACK_NUMERIC_REALLOC: " << Info[UMFPACK_NUMERIC_REALLOC] << "\n";
		std::cout << "    UMFPACK_NUMERIC_COSTLY_REALLOC: " << Info[UMFPACK_NUMERIC_COSTLY_REALLOC] << "\n";
		std::cout << "    UMFPACK_COMPRESSED_PATTERN: " << Info[UMFPACK_COMPRESSED_PATTERN] << "\n";
		std::cout << "    UMFPACK_LU_ENTRIES: " << Info[UMFPACK_LU_ENTRIES] << "\n";
		std::cout << "    UMFPACK_NUMERIC_TIME: " << Info[UMFPACK_NUMERIC_TIME] << "\n";
		std::cout << "    UMFPACK_RCOND: " << Info[UMFPACK_RCOND] << "\n";
		std::cout << "    UMFPACK_UDIAG_NZ: " << Info[UMFPACK_UDIAG_NZ] << "\n";
		std::cout << "    UMFPACK_UMIN: " << Info[UMFPACK_UMIN] << "\n";
		std::cout << "    UMFPACK_UMAX: " << Info[UMFPACK_UMAX] << "\n";
		std::cout << "    UMFPACK_MAX_FRONT_NROWS: " << Info[UMFPACK_MAX_FRONT_NROWS] << "\n";
		std::cout << "    UMFPACK_MAX_FRONT_NCOLS: " << Info[UMFPACK_MAX_FRONT_NCOLS] << "\n";
		std::cout << "    UMFPACK_ALL_LNZ: " << Info[UMFPACK_ALL_LNZ] << "\n";
		std::cout << "    UMFPACK_ALL_UNZ: " << Info[UMFPACK_ALL_UNZ] << "\n";
#endif
	    }

	    void init()
	    {
		MTL_THROW_IF(num_rows(A) != num_cols(A), matrix_not_square());
		n= int(num_rows(A));
		assign_pointers();
		init_aux(blong());
	    }



	  public:
	    /// Constructor referring to matrix \p A (not changed) and optionally Umfpack's strategy and alloc_init
	    solver(const matrix_type& A, int strategy = UMFPACK_STRATEGY_AUTO, double alloc_init = 0.7) 
	      : A(A), Apc(0), Aic(0), my_nnz(0), Symbolic(0), Numeric(0) 
	    {
		vampir_trace<5060> trace;
		// Use default setings.
		if (long_indices)
		    umfpack_dl_defaults(Control);
		else
		    umfpack_di_defaults(Control);

		Control[UMFPACK_STRATEGY] = strategy;
		Control[UMFPACK_ALLOC_INIT] = alloc_init;
		init(); 
	    }

	    ~solver()
	    {
		vampir_trace<5061> trace;
		if (long_indices) {
		    umfpack_dl_free_numeric(&Numeric);
		    umfpack_dl_free_symbolic(&Symbolic);
		} else {
		    umfpack_di_free_numeric(&Numeric);
		    umfpack_di_free_symbolic(&Symbolic);
		}
		if (Apc) delete[] Apc; 
		if (Aic) delete[] Aic;
	    }

	    void update_numeric_aux(true_)
	    {
		umfpack_dl_free_numeric(&Numeric);
		check(umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info), "Error in dl_numeric");
	    }
	    
	    void update_numeric_aux(false_)
	    {
		umfpack_di_free_numeric(&Numeric);
		check(umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info), "Error in di_numeric");
	    }

	    /// Update numeric part, for matrices that kept the sparsity and changed the values
	    void update_numeric()
	    {
		assign_pointers();
		update_numeric_aux(blong());
	    }

	    /// Update symbolic and numeric part
	    void update()
	    {
		if (long_indices) {
		    umfpack_dl_free_numeric(&Numeric);
		    umfpack_dl_free_symbolic(&Symbolic);
		} else {
		    umfpack_di_free_numeric(&Numeric);
		    umfpack_di_free_symbolic(&Symbolic);
		}
		init();
	    }

	    template <typename VectorX, typename VectorB>
	    void solve_aux(int sys, VectorX& xx, const VectorB& bb, true_)
	    {
	      check(umfpack_dl_solve(sys, Ap, Ai, Ax, &xx.value[0], &bb.value[0], Numeric, Control, Info), "Error in dl_solve");
	    }

	    template <typename VectorX, typename VectorB>
	    void solve_aux(int sys, VectorX& xx, const VectorB& bb, false_)
	    {
		check(umfpack_di_solve(sys, Ap, Ai, Ax, &xx.value[0], &bb.value[0], Numeric, Control, Info), "Error in di_solve");
	    }

	    /// Solve double system
	    template <typename VectorX, typename VectorB>
	    int operator()(VectorX& x, const VectorB& b)
	    {
		vampir_trace<5062> trace;
		MTL_THROW_IF(num_rows(A) != size(x) || num_rows(A) != size(b), incompatible_size());
		make_in_out_copy_or_reference<dense_vector<value_type>, VectorX> xx(x);
		make_in_copy_or_reference<dense_vector<value_type>, VectorB>     bb(b);
		int sys= mtl::traits::is_row_major<Parameters>::value ? UMFPACK_At : UMFPACK_A;
		solve_aux(sys, xx, bb, blong());
		return UMFPACK_OK;
	    }

	    /// Solve double system
	    template <typename VectorB, typename VectorX>
	    int solve(const VectorB& b, VectorX& x) const
	    {
		// return (*this)(x, b);
		return const_cast<solver&>(*this)(x, b); // evil hack because Umfpack has no const
	    }

	  private:
	    const matrix_type&  A;
	    int                 n;
	    const index_type    *Ap, *Ai;
	    index_type          *Apc, *Aic;
	    size_type           my_nnz;
	    const double        *Ax;
	    double              Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
	    void                *Symbolic, *Numeric;
	};

	/// Speciatization of solver for \ref mat::compressed2D with double values
	template <typename Parameters>
	class solver<compressed2D<std::complex<double>, Parameters> >
	  : public base_solver
	{
	    typedef std::complex<double>                      value_type;
	    typedef compressed2D<value_type, Parameters>      matrix_type;
	    typedef typename matrix_type::size_type           size_type;
	    typedef typename index<size_type>::type           index_type;

	    static const bool copy_indices= sizeof(index_type) != sizeof(size_type),
		              long_indices= use_long<size_type>::value;

	    typedef boost::mpl::bool_<long_indices>           blong;
	    typedef boost::mpl::true_                         true_;
	    typedef boost::mpl::false_                        false_;
	    

	    void assign_pointers()
	    {
		if (copy_indices) {
		    if (Apc == 0) Apc= new index_type[n + 1]; 
		    if (Aic == 0) Aic= new index_type[A.nnz()];
		    std::copy(A.address_major(), A.address_major() + n + 1, Apc);
		    std::copy(A.address_minor(), A.address_minor() + A.nnz(), Aic);
		    Ap= Apc;
		    Ai= Aic;
		} else {
		    Ap= reinterpret_cast<const index_type*>(A.address_major());
		    Ai= reinterpret_cast<const index_type*>(A.address_minor());
		}
		split_complex_vector(A.data, Ax, Az);
	    }

	    void init_aux(true_)
	    {
		check(umfpack_zl_symbolic(n, n, Ap, Ai, &Ax[0], &Az[0], &Symbolic, Control, Info), "Error in zl_symbolic");
		check(umfpack_zl_numeric(Ap, Ai, &Ax[0], &Az[0], Symbolic, &Numeric, Control, Info), "Error in zl_numeric");
	    }
	    
	    void init_aux(false_)
	    {
		check(umfpack_zi_symbolic(n, n, Ap, Ai, &Ax[0], &Az[0], &Symbolic, Control, Info), "Error in zi_symbolic");
		check(umfpack_zi_numeric(Ap, Ai, &Ax[0], &Az[0], Symbolic, &Numeric, Control, Info), "Error in zi_numeric");
	    }

	    void initialize()
	    {
		MTL_THROW_IF(num_rows(A) != num_cols(A), matrix_not_square());
		n= int(num_rows(A));
		assign_pointers();
		init_aux(blong());
	    }
	public:
	    /// Constructor referring to matrix \p A (not changed) and optionally Umfpack's strategy and alloc_init (look for the specializations)
	    explicit solver(const compressed2D<value_type, Parameters>& A, int strategy = UMFPACK_STRATEGY_AUTO, double alloc_init = 0.7) 
	      : A(A), Apc(0), Aic(0)
	    {
		vampir_trace<5060> trace;
		// Use default setings.
		if (long_indices)
		    umfpack_zl_defaults(Control);
		else
		    umfpack_zi_defaults(Control);
		// umfpack_zi_defaults(Control);

		Control[UMFPACK_STRATEGY] = strategy;
		Control[UMFPACK_ALLOC_INIT] = alloc_init;
		initialize();
	    }

	    ~solver()
	    {
		vampir_trace<5061> trace;
		if (long_indices) {
		    umfpack_zl_free_numeric(&Numeric);
		    umfpack_zl_free_symbolic(&Symbolic);
		} else {
		    umfpack_zi_free_numeric(&Numeric);
		    umfpack_zi_free_symbolic(&Symbolic);
		}
		if (Apc) delete[] Apc; 
		if (Aic) delete[] Aic;
	    }

	    void update_numeric_aux(true_)
	    {
		umfpack_zl_free_numeric(&Numeric);
		check(umfpack_zl_numeric(Ap, Ai, &Ax[0], &Az[0], Symbolic, &Numeric, Control, Info), "Error in dl_numeric D");
	    }
	    
	    void update_numeric_aux(false_)
	    {
		umfpack_zi_free_numeric(&Numeric);
		check(umfpack_zi_numeric(Ap, Ai, &Ax[0], &Az[0], Symbolic, &Numeric, Control, Info), "Error in di_numeric");
	    }

	    /// Update numeric part, for matrices that kept the sparsity and changed the values
	    void update_numeric()
	    {
		assign_pointers();
		update_numeric_aux(blong());
	    }

	    /// Update symbolic and numeric part
	    void update()
	    {
		Ax.change_dim(0); Az.change_dim(0);
		if (long_indices) {
		    umfpack_zl_free_numeric(&Numeric);
		    umfpack_zl_free_symbolic(&Symbolic);
		} else {
		    umfpack_zi_free_numeric(&Numeric);
		    umfpack_zi_free_symbolic(&Symbolic);
		}
		initialize();
	    }

	    template <typename VectorX, typename VectorB>
	    void solve_aux(int sys, VectorX& Xx, VectorX& Xz, const VectorB& Bx, const VectorB& Bz, true_)
	    {
		check(umfpack_zl_solve(sys, Ap, Ai, &Ax[0], &Az[0], &Xx[0], &Xz[0], &Bx[0], &Bz[0], Numeric, Control, Info), 
		      "Error in zi_solve");
	    }

	    template <typename VectorX, typename VectorB>
	    void solve_aux(int sys, VectorX& Xx, VectorX& Xz, const VectorB& Bx, const VectorB& Bz, false_)
	    {
		check(umfpack_zi_solve(sys, Ap, Ai, &Ax[0], &Az[0], &Xx[0], &Xz[0], &Bx[0], &Bz[0], Numeric, Control, Info), 
		      "Error in zi_solve");
	    }

	    /// Solve complex system
	    template <typename VectorX, typename VectorB>
	    int operator()(VectorX& x, const VectorB& b)
	    {
		vampir_trace<5062> trace;
		MTL_THROW_IF(num_rows(A) != size(x) || num_rows(A) != size(b), incompatible_size());
		dense_vector<double> Xx(size(x)), Xz(size(x)), Bx, Bz;
		split_complex_vector(b, Bx, Bz);
		int sys= mtl::traits::is_row_major<Parameters>::value ? UMFPACK_Aat : UMFPACK_A;
		solve_aux(sys, Xx, Xz, Bx, Bz, blong());
		merge_complex_vector(Xx, Xz, x);
		return UMFPACK_OK;
	    }

	    /// Solve complex system
	    template <typename VectorB, typename VectorX>
	    int solve(const VectorB& b, VectorX& x)
	    {
		return (*this)(x, b);
	    }

	private:
	    const matrix_type&   A;
	    int                  n; 
	    const index_type     *Ap, *Ai;
	    index_type          *Apc, *Aic;
	    dense_vector<double> Ax, Az;
	    double               Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
	    void                 *Symbolic, *Numeric;
	};

	template <typename Value, typename Parameters>
	class solver<compressed2D<Value, Parameters> >
	  : matrix_copy<compressed2D<Value, Parameters>, Value, typename Parameters::orientation>,
	    public solver<typename matrix_copy<compressed2D<Value, Parameters>, Value, typename Parameters::orientation>::matrix_type >
	{
	    typedef matrix_copy<compressed2D<Value, Parameters>, Value, typename Parameters::orientation> copy_type;
	    typedef solver<typename matrix_copy<compressed2D<Value, Parameters>, Value, typename Parameters::orientation>::matrix_type > solver_type;
	public:
	    explicit solver(const compressed2D<Value, Parameters>& A) 
		: copy_type(A), solver_type(copy_type::matrix), A(A)
	    {}

	    void update()
	    {
		copy_type::matrix= A;
		solver_type::update();
	    }

	    void update_numeric()
	    {
		copy_type::matrix= A;
		solver_type::update_numeric();
	    }
	private:
	    const compressed2D<Value, Parameters>& A;
	};


    } // umfpack

/// Convenience function to create a solver
template <typename Value, typename Parameters>
umfpack::solver<compressed2D<Value, Parameters> > 
inline make_umfpack_solver(const compressed2D<Value, Parameters>& A)
{  return umfpack::solver<compressed2D<Value, Parameters> >(A); }


/// Solve A*x == b with umfpack
/** Only available when compiled with enabled macro MTL_HAS_UMFPACK.
    Uses classes umfpack::solver internally.
    If you want more control on single operations or to keep umfpack's
    internal factorization, use this class.
 **/
template <typename Value, typename Parameters, typename VectorX, typename VectorB>
int umfpack_solve(const compressed2D<Value, Parameters>& A, VectorX& x, const VectorB& b)
{
    umfpack::solver<compressed2D<Value, Parameters> > solver(A);
    return solver(x, b);
}

}} // namespace mtl::mat

#endif

#endif // MTL_MATRIX_UMFPACK_SOLVE_INCLUDE
