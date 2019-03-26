

extern "C" {
#  include <umfpack.h>
}


#if 0 // for checking only presence of libraries to link not includes
extern "C" {
int umfpack_di_symbolic
     (int n_row,
      int n_col,
      const int Ap [ ],
      const int Ai [ ],
      const double Ax [ ],
      void **Symbolic,
      const double Control [3],
      double Info [4]
      );
}
#endif

int main()
{
    int    ai[]= {1, 2};
    double ad[]= {1.0, 2.0};
    void*  vp;
    
    umfpack_di_symbolic(1, 2, ai, ai, ad, &vp, ad, ad);
}
