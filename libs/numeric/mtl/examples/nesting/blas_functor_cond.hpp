template <typename MatrixA, typename MatrixB, typename MatrixC>
struct blas_mult_ft
    : public mult_ft<MatrixA, MatrixB, MatrixC>
{};

#ifdef MTL_HAS_BLAS

/* ... here come the specializations */

#endif // MTL_HAS_BLAS
