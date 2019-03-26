using assign::assign_sum;
typedef gen_platform_dmat_dmat_mult_t<assign_sum, gen_tiling_44_dmat_dmat_mult_t> platform_mult_type;
gen_blas_dmat_dmat_mult_t<assign_sum, platform_mult_type>                         my_mult;

// ...
my_mult(A, B, C);
