template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_blas_dmat_dmat_mult_ft
    : public Backup
{};

/* ... its specializations */

template <typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_blas_dmat_dmat_mult_t
    : public Backup
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& a, MatrixB const& b, MatrixC& c)
    {
	gen_blas_dmat_dmat_mult_ft<MatrixA, MatrixB, MatrixC, Assign, Backup>()(a, b, c);
    }
};
