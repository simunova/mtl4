
/*
 *  baseCasesBoost.h
 *  Jan07 
 *  Created by  Adwait Joshi
 *  Project: Arcee
 *
 *  Copyright (c) 2007 The Trustees of Indiana University. All rights reserved.
 *
 */

 
long unsigned int docholeskyhits = 0;
long unsigned int schurhits = 0;
long unsigned int trischurhits = 0;
long unsigned int trisolvehits = 0;

template <typename Matrix>
void do_cholesky_base(Matrix& matrix)
{ 
  for(int k=0; k<matrix.num_rows();  k++ ){
     matrix[k][k] = sqrt( matrix[k][k]);     
      
     for(int i=k+1; i<matrix.num_rows(); i++ )
			 matrix[i][k] /= matrix[k][k];
				
     for(int i=k+1; i<matrix.num_rows();  i++){
     	 typename Matrix::value_type d = matrix[i][k];
			 for (int j=k+1; j<=i; j++ ) 
				  matrix[i][j] -= d * matrix[j][k];		    
	  }
	}
  return;
}


template <typename Matrix>
void tri_solve_base(Matrix& SW ,Matrix& NW)
{
	for(int k=0; k<NW.num_rows();  k++ ){

		for(int i=0; i<SW.num_rows(); i++ )
				SW[i][k] /= NW[k][k];
		
		for(int i=0 ; i<SW.num_rows();  i++ ){
			 typename Matrix::value_type d = SW[i][k];
			 for (int j=k+1; j<SW.num_cols(); j++ ) 
				SW[i][j] -= d * NW[j][k];			 
		}
	}
	return;
}


template <typename MatrixSE, typename MatrixSW>
void tri_schur_base(MatrixSE& SE ,MatrixSW& SW)
{		
 	for(int k=0; k<SW.num_rows();  k++ ) 
		for (int i=0 ; i < SE.num_rows();  i++ )	{
			typename MatrixSW::value_type d = SW[i][k];
			for (int j=0; j<=i; j++ ) 
				SE[i][j] -= d * SW[j][k];
		}
}




#if 0
template <typename MatrixSE, typename MatrixSW>
struct tri_schur_base_bracket_t
{
	void operator() (MatrixSE& SE ,MatrixSW& SW)
	{		
		//trischurhits++;
		for(int k=0; k<SW.num_rows();  k++ ) {

			for (int i=0 ; i < SE.num_rows();  i++ )	{
				typename MatrixSW::value_type d = SW[i][k];
				for (int j=0; j<=i; j++ ) {
 					SE[i][j] -= d * SW[j][k];
 				}
			}
		}
	}
};

template <typename TriSchurBase, ...> 
struct rec_cholesky { .... 
	
	TriSchurBase  tsb; ... tsb(se, sw);
	TriSchurBase()(se, sw);
	}


rec_cholesky<tri_schur_base_bt<Matrix, Matrix>, ...>

schur_base(_1, _2, _3)

#endif





template <typename Matrix>
void schur_base(Matrix& NE ,Matrix& NW ,Matrix& SW)
{ 
	for(int k=0; k<NW.num_rows();  k++ ){
		for(int i=0 ; i<NE.num_rows();  i++ ){
			typename Matrix::value_type d = NW[i][k];
			for (int j=0; j<NE.num_cols(); j++ )			 
				NE[i][j] -= d * SW[j][k];				 
		}
	}
	return;
}



template <typename Matrix>
void verify_matrix(Matrix& A)
{ 
	
	
	
	for(int i=0; i<A.num_rows();  i++ ){
		for(int j=0 ; j<A.num_rows();  j++ ){
		
			typename Matrix::value_type  element=0;
			for (int k=0; k<A.num_cols(); k++ )			 
				element += A[i][k] * A[j][k];		
		cout<<element<<"  ";		 
		}
		cout<<endl;
	}
	return;
}


 
