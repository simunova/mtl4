template <typename Backup>
struct blas_mult_t
    : public Backup
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& a, MatrixB const& b, MatrixC& c) 
    { 
	blas_mult_ft<MatrixA, MatrixB, MatrixC, Backup>()(a, b, c);
    }
};
