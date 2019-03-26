#include <boost/numeric/mtl/interface/vpt.hpp>

#ifdef MTL_HAS_VPT
namespace mtl { 

/// Namespace for Vampir Trace interface
namespace vpt {

// Categories:
// Utilities + very small functions:  0000
// Static size operations:            1000
// Vector operations:                 2000
// Matrix Vector & single matrix:     3000
// Matrix matrix operations:          4000
// Factorizations, preconditioners:   5000
// Fused operations:                  6000
// Iterative solvers:                 7000
// Multigrid:                         8000


// Utilities:                        < 1000
template <> std::string vampir_trace<1>::name("copysign");
template <> std::string vampir_trace<2>::name("Elem_raw_copy");
template <> std::string vampir_trace<3>::name("Get_real_part");
template <> std::string vampir_trace<4>::name("Info_contruct_vector");
template <> std::string vampir_trace<5>::name("right_scale_inplace");
template <> std::string vampir_trace<6>::name("sign_real_part_of_complex");
template <> std::string vampir_trace<7>::name("unrolling_expresion");
template <> std::string vampir_trace<8>::name("");
template <> std::string vampir_trace<9>::name("");
template <> std::string vampir_trace<10>::name("squared_abs_magnitudes");
template <> std::string vampir_trace<11>::name("squared_abs_complex");
template <> std::string vampir_trace<12>::name("squared_abs_magnitudes_template");
template <> std::string vampir_trace<13>::name("update_store");
template <> std::string vampir_trace<14>::name("update_plus");
template <> std::string vampir_trace<15>::name("update_minus");
template <> std::string vampir_trace<16>::name("update_times");
template <> std::string vampir_trace<17>::name("update_adapter");
template <> std::string vampir_trace<18>::name("");
template <> std::string vampir_trace<19>::name("");
template <> std::string vampir_trace<20>::name("update_proxy_<<");
template <> std::string vampir_trace<21>::name("update_proxy_=");
template <> std::string vampir_trace<22>::name("update_proxy_+=");
template <> std::string vampir_trace<23>::name("sfunctor::plus");
template <> std::string vampir_trace<24>::name("sfunctor::minus");
template <> std::string vampir_trace<25>::name("sfunctor::times");
template <> std::string vampir_trace<26>::name("sfunctor::divide");
template <> std::string vampir_trace<27>::name("sfunctor::assign");
template <> std::string vampir_trace<28>::name("sfunctor::plus_assign");
template <> std::string vampir_trace<29>::name("sfunctor::minus_assign");
template <> std::string vampir_trace<30>::name("sfunctor::times_assign");
template <> std::string vampir_trace<31>::name("sfunctor::divide_assign");
template <> std::string vampir_trace<32>::name("sfunctor::identity");
template <> std::string vampir_trace<33>::name("sfunctor::abs");
template <> std::string vampir_trace<34>::name("sfunctor::sqrt");
template <> std::string vampir_trace<35>::name("sfunctor::square");
template <> std::string vampir_trace<36>::name("sfunctor::negate");
template <> std::string vampir_trace<37>::name("sfunctor::compose");
template <> std::string vampir_trace<38>::name("sfunctor::compose_first");
template <> std::string vampir_trace<39>::name("sfunctor::compose_second");
template <> std::string vampir_trace<40>::name("sfunctor::compose_both");
template <> std::string vampir_trace<41>::name("sfunctor::compose_binary");
template <> std::string vampir_trace<42>::name("");
template <> std::string vampir_trace<43>::name("");
template <> std::string vampir_trace<44>::name("");
template <> std::string vampir_trace<45>::name("");

// Fine-grained vector operations
template <> std::string vampir_trace<236>::name("Vector_swapped_row");


// Static size operations:           1000
template <> std::string vampir_trace<1001>::name("stat_vec_expr");
template <> std::string vampir_trace<1002>::name("fsize_dmat_dmat_mult");
template <> std::string vampir_trace<1003>::name("vector_size_static");
template <> std::string vampir_trace<1004>::name("static_dispatch"); // ?? row_in_matrix.hpp:74
template <> std::string vampir_trace<1005>::name("copy_blocks_forward");
template <> std::string vampir_trace<1006>::name("copy_blocks_backward");
template <> std::string vampir_trace<1007>::name("Static_Size");
template <> std::string vampir_trace<1008>::name("fsize_mat_vect_mult");
template <> std::string vampir_trace<1009>::name("");
template <> std::string vampir_trace<1010>::name("");
template <> std::string vampir_trace<1011>::name("");
template <> std::string vampir_trace<1012>::name("");
template <> std::string vampir_trace<1013>::name("");
template <> std::string vampir_trace<1014>::name("");
template <> std::string vampir_trace<1015>::name("");
template <> std::string vampir_trace<1016>::name("");
template <> std::string vampir_trace<1017>::name("");
template <> std::string vampir_trace<1018>::name("");
template <> std::string vampir_trace<1019>::name("");
template <> std::string vampir_trace<1020>::name("");





// Vector operations:                2000
template <> std::string vampir_trace<2001>::name("gen_vector_copy");
template <> std::string vampir_trace<2002>::name("cross");
template <> std::string vampir_trace<2003>::name("dot");
template <> std::string vampir_trace<2004>::name("householder");
template <> std::string vampir_trace<2005>::name("householder_s");
template <> std::string vampir_trace<2006>::name("infinity_norm");
template <> std::string vampir_trace<2007>::name("look_at_each_nonzero");
template <> std::string vampir_trace<2008>::name("look_at_each_nonzero_pos");
template <> std::string vampir_trace<2009>::name("reduction");
template <> std::string vampir_trace<2010>::name("max");
template <> std::string vampir_trace<2011>::name("max_abs_pos");
template <> std::string vampir_trace<2012>::name("max_of_sums");
template <> std::string vampir_trace<2013>::name("max_pos");
template <> std::string vampir_trace<2014>::name("merge_complex_vector");
template <> std::string vampir_trace<2015>::name("one_norm");
template <> std::string vampir_trace<2016>::name("diagonal");
template <> std::string vampir_trace<2017>::name("dyn_vec_expr");
template <> std::string vampir_trace<2018>::name("Orthogonalize_Vectors");
template <> std::string vampir_trace<2019>::name("Orthogonalize_Factors");				      
template <> std::string vampir_trace<2020>::name("Vector_product");
template <> std::string vampir_trace<2021>::name("Vector_random");
template <> std::string vampir_trace<2022>::name("Vec_Vec_rank_update");
template <> std::string vampir_trace<2023>::name("Vector_dispatch");
template <> std::string vampir_trace<2024>::name("Vector_rscale");
template <> std::string vampir_trace<2025>::name("Multi-vector_mult");
template <> std::string vampir_trace<2026>::name("Transp_Multi-vector_mult");
template <> std::string vampir_trace<2027>::name("Hermitian_Multi-vector_mult");
template <> std::string vampir_trace<2028>::name("Vector_scal");
template <> std::string vampir_trace<2029>::name("Vector_set_zero");
template <> std::string vampir_trace<2030>::name("Vector_size1D");
template <> std::string vampir_trace<2031>::name("Vector_size_runtime");
template <> std::string vampir_trace<2032>::name("Vect_quicksort_lo_to_hi");
template <> std::string vampir_trace<2033>::name("Vect_quicksort_permutaion_lo_to_hi");
template <> std::string vampir_trace<2034>::name("split_complex_vector");
template <> std::string vampir_trace<2035>::name("Vect_entries_sum");
template <> std::string vampir_trace<2037>::name("Vector_const_trans");
template <> std::string vampir_trace<2038>::name("Vector_trans");
template <> std::string vampir_trace<2039>::name("two_norm");
template <> std::string vampir_trace<2040>::name("dot_simple");
template <> std::string vampir_trace<2041>::name("unary_dot");
template <> std::string vampir_trace<2042>::name("dense_copy_ctor");
template <> std::string vampir_trace<2043>::name("dense_tpl_copy_ctor");
template <> std::string vampir_trace<2044>::name("");
template <> std::string vampir_trace<2045>::name("");
template <> std::string vampir_trace<2046>::name("");
template <> std::string vampir_trace<2047>::name("");
template <> std::string vampir_trace<2048>::name("");
template <> std::string vampir_trace<2049>::name("");
template <> std::string vampir_trace<2050>::name("");
template <> std::string vampir_trace<2051>::name("");
template <> std::string vampir_trace<2052>::name("");


// Matrix Vector & single matrix:    3000
template <> std::string vampir_trace<3001>::name("matrix_copy_ele_times");
template <> std::string vampir_trace<3002>::name("gen_matrix_copy");
template <> std::string vampir_trace<3003>::name("copy");
template <> std::string vampir_trace<3004>::name("clone");
template <> std::string vampir_trace<3005>::name("compute_summand");
template <> std::string vampir_trace<3006>::name("crop");
template <> std::string vampir_trace<3007>::name("mat::diagonal");
template <> std::string vampir_trace<3008>::name("assign_each_nonzero");
template <> std::string vampir_trace<3009>::name("fill");
template <> std::string vampir_trace<3010>::name("frobenius_norm");
template <> std::string vampir_trace<3011>::name("mat::infinity_norm");
template <> std::string vampir_trace<3012>::name("invert_diagonal");
template <> std::string vampir_trace<3013>::name("iota");
template <> std::string vampir_trace<3014>::name("left_scale_inplace");
template <> std::string vampir_trace<3015>::name("mat::look_at_each_nonzero");
template <> std::string vampir_trace<3016>::name("mat::look_at_each_nonzero_pos");
template <> std::string vampir_trace<3017>::name("fsize_dense_mat_cvec_mult");
template <> std::string vampir_trace<3018>::name("dense_mat_cvec_mult");
template <> std::string vampir_trace<3019>::name("mvec_cvec_mult");
template <> std::string vampir_trace<3020>::name("trans_mvec_cvec_mult");
template <> std::string vampir_trace<3021>::name("herm_mvec_cvec_mult");
template <> std::string vampir_trace<3022>::name("sparse_row_cvec_mult"); // generic row-major sparse
template <> std::string vampir_trace<3023>::name("ccs_cvec_mult");
template <> std::string vampir_trace<3024>::name("mat::max_abs_pos");
template <> std::string vampir_trace<3025>::name("mat::one_norm");
template <> std::string vampir_trace<3026>::name("invert_diagonal(compressed)");
template <> std::string vampir_trace<3027>::name("mat_vect_mult");
template <> std::string vampir_trace<3028>::name("Vect_sparse_mat_mult");
template <> std::string vampir_trace<3029>::name("Matrix_scal");
template <> std::string vampir_trace<3030>::name("Vector_Secular_Equation");
template <> std::string vampir_trace<3031>::name("Matrix_set_zero");
template <> std::string vampir_trace<3032>::name("Matrix_size1D");
template <> std::string vampir_trace<3033>::name("Matrix_size_runtime");
template <> std::string vampir_trace<3034>::name("Matrix_LU");
template <> std::string vampir_trace<3035>::name("Vector_Matrix_LU");
template <> std::string vampir_trace<3036>::name("Sub_matrix_indices");
template <> std::string vampir_trace<3037>::name("Matrix_svd_reference");
template <> std::string vampir_trace<3038>::name("Matrix_svd_triplet");
template <> std::string vampir_trace<3039>::name("Matrix_swapped");
template <> std::string vampir_trace<3040>::name("Matrix_Trace");
template <> std::string vampir_trace<3041>::name("Matrix_const_trans");
template <> std::string vampir_trace<3042>::name("Matrix_trans");
template <> std::string vampir_trace<3043>::name("Matrix_upper_trisolve");
template <> std::string vampir_trace<3044>::name("Matrix_upper_trisolve_diagonal");
template <> std::string vampir_trace<3045>::name("Matrix_upper_trisolve_invers_diag");
template <> std::string vampir_trace<3046>::name("Matrix_upper_trisolve_DiaTag");
template <> std::string vampir_trace<3047>::name("scalar_assign");
template <> std::string vampir_trace<3048>::name("elest_cvec_mult");
template <> std::string vampir_trace<3049>::name("crs_cvec_mult");
template <> std::string vampir_trace<3050>::name("sparse_ins::ctor");
template <> std::string vampir_trace<3051>::name("sparse_ins::dtor");
template <> std::string vampir_trace<3052>::name("sparse_ins::stretch");
template <> std::string vampir_trace<3053>::name("sparse_ins::final_place");
template <> std::string vampir_trace<3054>::name("sparse_ins::insert_spare");
template <> std::string vampir_trace<3055>::name("mat_crtp_scal_assign");
template <> std::string vampir_trace<3056>::name("mat_crtp_mat_assign");
template <> std::string vampir_trace<3057>::name("mat_crtp_sum_assign");
template <> std::string vampir_trace<3058>::name("mat_crtp_diff_assign");
template <> std::string vampir_trace<3059>::name("mat_crtp_array_assign");
template <> std::string vampir_trace<3060>::name("mat_crtp_mvec_assign");
template <> std::string vampir_trace<3061>::name("copy_band_to_sparse");
template <> std::string vampir_trace<3062>::name("block_dia_times_cvec");
template <> std::string vampir_trace<3063>::name("laplacian_setup");
template <> std::string vampir_trace<3064>::name("vsmat_cvec_mult");
template <> std::string vampir_trace<3065>::name("adapt_crs_cvec_mult");
template <> std::string vampir_trace<3066>::name("dense2D_cvec_mult");
template <> std::string vampir_trace<3067>::name("square_cvec_mult");
template <> std::string vampir_trace<3068>::name("mat_crtp_mult_assign");
template <> std::string vampir_trace<3069>::name("sbanded_cvec_mult");
template <> std::string vampir_trace<3070>::name("mat_cvec_multiplier");
template <> std::string vampir_trace<3071>::name("");
template <> std::string vampir_trace<3072>::name("");
template <> std::string vampir_trace<3073>::name("");
template <> std::string vampir_trace<3074>::name("");
template <> std::string vampir_trace<3075>::name("");
template <> std::string vampir_trace<3076>::name("");
template <> std::string vampir_trace<3077>::name("");
template <> std::string vampir_trace<3078>::name("");
template <> std::string vampir_trace<3079>::name("");


// Matrix matrix operations:        4000
template <> std::string vampir_trace<4001>::name("cursor_dmat_dmat_mult");
template <> std::string vampir_trace<4002>::name("dmat_dmat_mult");
template <> std::string vampir_trace<4003>::name("tiling_dmat_dmat_mult");
template <> std::string vampir_trace<4004>::name("tiling_44_dmat_dmat_mult");
template <> std::string vampir_trace<4005>::name("tiling_22_dmat_dmat_mult");
template <> std::string vampir_trace<4006>::name("wrec_dmat_dmat_mult");
template <> std::string vampir_trace<4007>::name("recursive_dmat_dmat_mult");
template <> std::string vampir_trace<4008>::name("xgemm");
template <> std::string vampir_trace<4009>::name("");
template <> std::string vampir_trace<4010>::name("mult");
template <> std::string vampir_trace<4011>::name("gen_mult");
template <> std::string vampir_trace<4012>::name("mat_mat_mult");
template <> std::string vampir_trace<4013>::name("matrix_qr");
template <> std::string vampir_trace<4014>::name("matrix_qr_factors");
template <> std::string vampir_trace<4015>::name("matrix_random");
template <> std::string vampir_trace<4016>::name("matrix_scale_inplace");
template <> std::string vampir_trace<4017>::name("matrix_rscale");
template <> std::string vampir_trace<4018>::name("matrix_gen_smat_dmat_mult");
template <> std::string vampir_trace<4019>::name("matrix_gen_tiling_smat_dmat_mult");    
template <> std::string vampir_trace<4020>::name("matrix_smat_smat_mult");
template <> std::string vampir_trace<4021>::name("");
template <> std::string vampir_trace<4022>::name("");
template <> std::string vampir_trace<4023>::name("");
template <> std::string vampir_trace<4024>::name("");
template <> std::string vampir_trace<4025>::name("");
template <> std::string vampir_trace<4026>::name("");
template <> std::string vampir_trace<4027>::name("");
template <> std::string vampir_trace<4028>::name("");
template <> std::string vampir_trace<4029>::name("");
template <> std::string vampir_trace<4030>::name("");
template <> std::string vampir_trace<4031>::name("");
template <> std::string vampir_trace<4032>::name("");
template <> std::string vampir_trace<4033>::name("");
template <> std::string vampir_trace<4034>::name("");
template <> std::string vampir_trace<4035>::name("");
template <> std::string vampir_trace<4036>::name("read_el_matrix");
template <> std::string vampir_trace<4037>::name("");
template <> std::string vampir_trace<4038>::name("");
template <> std::string vampir_trace<4039>::name("");
template <> std::string vampir_trace<4040>::name("");
template <> std::string vampir_trace<4041>::name("");


// Factorizations, preconditioners: 5000
template <> std::string vampir_trace<5001>::name("cholesky_base");
template <> std::string vampir_trace<5002>::name("cholesky_solve_base");
template <> std::string vampir_trace<5003>::name("cholesky_schur_base");
template <> std::string vampir_trace<5004>::name("cholesky_update_base");
template <> std::string vampir_trace<5005>::name("cholesky_schur_update");
template <> std::string vampir_trace<5006>::name("cholesky_tri_solve");
template <> std::string vampir_trace<5007>::name("cholesky_tri_schur");
template <> std::string vampir_trace<5008>::name("recursive cholesky");
template <> std::string vampir_trace<5009>::name("fill_matrix_for_cholesky");
template <> std::string vampir_trace<5010>::name("qr_sym_imp");
template <> std::string vampir_trace<5011>::name("qr_algo");
template <> std::string vampir_trace<5012>::name("eigenvalue_symmetric");
template <> std::string vampir_trace<5013>::name("hessenberg_q");
template <> std::string vampir_trace<5014>::name("hessenberg_factors");
template <> std::string vampir_trace<5015>::name("extract_householder_hessenberg");
template <> std::string vampir_trace<5016>::name("householder_hessenberg");
template <> std::string vampir_trace<5017>::name("extract_hessenberg");
template <> std::string vampir_trace<5018>::name("hessenberg");
template <> std::string vampir_trace<5019>::name("inv_upper");
template <> std::string vampir_trace<5020>::name("inv_lower");
template <> std::string vampir_trace<5021>::name("inv");
template <> std::string vampir_trace<5022>::name("lower_trisolve");
template <> std::string vampir_trace<5023>::name("lu");
template <> std::string vampir_trace<5024>::name("lu(pivot)");
template <> std::string vampir_trace<5025>::name("lu_f");
template <> std::string vampir_trace<5026>::name("lu_solve_straight");
template <> std::string vampir_trace<5027>::name("lu_apply");
template <> std::string vampir_trace<5028>::name("lu_solve");
template <> std::string vampir_trace<5029>::name("lu_adjoint_apply");
template <> std::string vampir_trace<5030>::name("lu_adjoint_solve");
template <> std::string vampir_trace<5031>::name("pc::id::solve");
template <> std::string vampir_trace<5032>::name("pc::id.solve");
template <> std::string vampir_trace<5033>::name("pc::id::adjoint_solve");
template <> std::string vampir_trace<5034>::name("pc::id.adjoint_solve");
template <> std::string vampir_trace<5035>::name("ic_0::factorize");
template <> std::string vampir_trace<5036>::name("ic_0::solve");
template <> std::string vampir_trace<5037>::name("ic_0::solve_nocopy");
template <> std::string vampir_trace<5038>::name("ilu_0::factorize");
template <> std::string vampir_trace<5039>::name("ilu_0::solve");
template <> std::string vampir_trace<5040>::name("ilu_0::adjoint_solve");
template <> std::string vampir_trace<5041>::name("lower_trisolve_kernel");
template <> std::string vampir_trace<5042>::name("upper_trisolve_row");
template <> std::string vampir_trace<5043>::name("upper_trisolve_col");
template <> std::string vampir_trace<5044>::name("ic_0::adjoint_solve");
template <> std::string vampir_trace<5045>::name("ic_0::adjoint_solve_nocopy");
template <> std::string vampir_trace<5046>::name("upper_trisolve_crs_compact");
template <> std::string vampir_trace<5047>::name("lower_trisolve_crs_compact");
template <> std::string vampir_trace<5048>::name("lower_unit_trisolve_crs_compact");
template <> std::string vampir_trace<5049>::name("ilut::factorize");
template <> std::string vampir_trace<5050>::name("diagonal::setup");
template <> std::string vampir_trace<5051>::name("diagonal::solve");
template <> std::string vampir_trace<5052>::name("imf::factor");
template <> std::string vampir_trace<5053>::name("imf::ctor");
template <> std::string vampir_trace<5054>::name("imf::solve");
template <> std::string vampir_trace<5055>::name("pc::solver::assign_to");
template <> std::string vampir_trace<5056>::name("sub_matrix_pc::solve");
template <> std::string vampir_trace<5057>::name("sub_matrix_pc::adjoint_solve");
template <> std::string vampir_trace<5058>::name("pc::concat::solve");
template <> std::string vampir_trace<5059>::name("pc::concat::adjoint_solve");
template <> std::string vampir_trace<5060>::name("umfpack::solver::ctor");
template <> std::string vampir_trace<5061>::name("umfpack::solver::dtor");
template <> std::string vampir_trace<5062>::name("umfpack::solve");


// Fused operations:                6000
template <> std::string vampir_trace<6001>::name("fused::fwd_eval_loop");
template <> std::string vampir_trace<6002>::name("fused::fwd_eval_loop_unrolled");
template <> std::string vampir_trace<6003>::name("fused::bwd_eval_loop");
template <> std::string vampir_trace<6004>::name("fused::bwd_eval_loop_unrolled");



// Iterative solvers:               7000
template <> std::string vampir_trace<7001>::name("cg_without_pc");
template <> std::string vampir_trace<7002>::name("cg");
template <> std::string vampir_trace<7003>::name("bicg");
template <> std::string vampir_trace<7004>::name("bicgstab");
template <> std::string vampir_trace<7005>::name("bicgstab_2");
template <> std::string vampir_trace<7006>::name("bicgstab_ell");
template <> std::string vampir_trace<7007>::name("cgs");
template <> std::string vampir_trace<7008>::name("qmr");
template <> std::string vampir_trace<7009>::name("tfqmr");
template <> std::string vampir_trace<7010>::name("idr_s");


// OpenMP
template <> std::string vampir_trace<8001>::name("omp::dot");
template <> std::string vampir_trace<8002>::name("omp::reduction");
template <> std::string vampir_trace<8003>::name("omp::dyn_vec_expr");
template <> std::string vampir_trace<8004>::name("omp::crs_cvec_mult");



// multigrid
template <> std::string vampir_trace<8501>::name("mtl::mg::v_cycle");
template <> std::string vampir_trace<8502>::name("mtl::mg::w_cycle");
template <> std::string vampir_trace<8503>::name("mtl::mg::fmg");
template <> std::string vampir_trace<8504>::name("mtl::mg::two_grid_cycle");

template <> std::string vampir_trace<8510>::name("mtl::mg::geometric_multigrid_solver_impl");
template <> std::string vampir_trace<8511>::name("mtl::mg::geometric_multigrid_solver_solve1");
template <> std::string vampir_trace<8512>::name("mtl::mg::geometric_multigrid_solver_solve2");

template <> std::string vampir_trace<8515>::name("mtl::mg::algebraic_multigrid_solver");
template <> std::string vampir_trace<8516>::name("amg_pc::solve");

template <> std::string vampir_trace<8520>::name("mtl::mg::linear_restriction");
template <> std::string vampir_trace<8521>::name("mtl::mg::linear_prolongation");

template <> std::string vampir_trace<8530>::name("mtl::mg::gauss_elimination");
template <> std::string vampir_trace<8531>::name("mtl::mg::back_substitution");

template <> std::string vampir_trace<8550>::name("mtl::mg::jacobi");
template <> std::string vampir_trace<8551>::name("mtl::mg::gauss_seidel");
template <> std::string vampir_trace<8552>::name("mtl::mg::jor");
template <> std::string vampir_trace<8553>::name("mtl::mg::sor");

template <> std::string vampir_trace<8572>::name("boundaries");
template <> std::string vampir_trace<8573>::name("viscosity");
template <> std::string vampir_trace<8574>::name("pressure_correction");

template <> std::string vampir_trace<8590>::name("mtl::mg::util::vtk_exporter");
template <> std::string vampir_trace<8591>::name("mtl::mg::util::csv_exporter");

template <> std::string vampir_trace<8610>::name("amg::amg_matrix_hierarchy");
template <> std::string vampir_trace<8611>::name("amg::compute_influence");
template <> std::string vampir_trace<8612>::name("amg::default_coarse_grid_detection::compute_C");
template <> std::string vampir_trace<8614>::name("amg::utils::compute_potentials");
template <> std::string vampir_trace<8615>::name("amg::utils::find_max_pos");

template <> std::string vampir_trace<8617>::name("amg::amg_prolongation");
template <> std::string vampir_trace<8618>::name("amg::compute_weight");
template <> std::string vampir_trace<8619>::name("amg::compute_mfactors");

template <> std::string vampir_trace<8620>::name("amg::strongly_influenced_points");
template <> std::string vampir_trace<8621>::name("amg::is_strongly_influenced");
template <> std::string vampir_trace<8622>::name("amg::strongly_influencing_points");

template <> std::string vampir_trace<8630>::name("amg::amg_operators::amg_restriction");
template <> std::string vampir_trace<8631>::name("amg::amg_operators::amg_prolongation");
template <> std::string vampir_trace<8635>::name("amg::amg_operators::amg_weight");

template <> std::string vampir_trace<8900>::name("NaSto::solve()");
template <> std::string vampir_trace<8910>::name("NaSto::computeGamma()");
template <> std::string vampir_trace<8920>::name("NaSto::computeBoundaries()");
template <> std::string vampir_trace<8930>::name("NaSto::computeImplViscosity()");
template <> std::string vampir_trace<8940>::name("NaSto::computePressureCorr()");

// Test blocks for performance debugging
template <> std::string vampir_trace<9901>::name("tb1"); 
template <> std::string vampir_trace<9902>::name("tb2");
template <> std::string vampir_trace<9903>::name("tb3");
template <> std::string vampir_trace<9904>::name("tb4");
template <> std::string vampir_trace<9999>::name("main");
    

// Only for testing
template <> std::string vampir_trace<9990>::name("helper_function");
template <> std::string vampir_trace<9991>::name("function");



}} //mtl::vpt

#endif //MTL_HAS_VPT
