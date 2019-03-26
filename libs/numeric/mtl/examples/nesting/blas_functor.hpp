template <typename MatrixA, typename MatrixB, typename MatrixC>
struct blas_mult_ft
    : public mult_ft<MatrixA, MatrixB, MatrixC>
{};

/* ... here come the specializations */

