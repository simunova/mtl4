template <typename MatrixA, typename MatrixB, typename MatrixC, 
	  typename Backup= mult_ft<MatrixA, MatrixB, MatrixC> >
struct blas_mult_ft
    : public Backup
{};
