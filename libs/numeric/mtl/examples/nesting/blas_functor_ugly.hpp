template <typename MatrixA, typename MatrixB, typename MatrixC>
struct blas_mult_ft
{
    void operator()(MatrixA const& a, MatrixB const& b, MatrixC& c) 
    { 
	mult_ft<MatrixA, MatrixB, MatrixC>()(a, b, c);
    }
};


template <typename ParaA, typename ParaB, typename ParaC>
struct blas_mult_ft<dense2D<double, ParaA>, dense2D<double, ParaB>, 
		    dense2D<double, ParaC> >
{
    void operator()(const dense2D<double, ParaA>& a, const dense2D<double, ParaB>& b, 
		    dense2D<double, ParaC>& c)
    {
	/* ... _dgemm( ... only 13 arguments ...); */
    }
};

/* ... more specializations */
