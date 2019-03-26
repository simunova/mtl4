/*
 *  strassen_rectWise.cpp
 *  July 17 
 * 
 *  Edited by Sarah Loos
 *  Project: Arcee
 *  
 *  Copyright (c) 2007-2008 The Trustees of Indiana University. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <boost/test/minimal.hpp>
#include <boost/tuple/tuple.hpp>
//#include <papi.h>
#include <vector>

#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/operation/sub_matrix.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>
// Replaced preceding file: recurator -> recursator.
#include <boost/numeric/mtl/recursion/base_case_test.hpp>
#include <boost/numeric/mtl/recursion/for_each.hpp>
#include <boost/numeric/mtl/recursion/base_case_cast.hpp>
#include <boost/numeric/mtl/recursion/base_case_matrix.hpp>
#include <boost/numeric/mtl/recursion/bit_masking.hpp>


using namespace mtl;
using namespace std; 
using namespace mtl::recursion; 

#define NUM_EVENTS 4

const unsigned long mask = generate_mask<true,  2 ,row_major, 0>::value;
typedef morton_dense<double,  mask> matrix_type;

// 'order' Can be passed as command line parameter, see main -- pg

const int basecasesize =4;
typedef recursion::bound_test_static<basecasesize> BaseTest;
BaseTest    is_base;

int  callnum = 0, basehit = 0;
int docholcall=0;
int schurcall=0;
int trischurcall=0;
int trisolvecall=0;


#define lgMtx 4

#define allocTempMtx(lev)                               \
matrix_type M0_##lev(1<<(lgMtx-lev-1),1<<(lgMtx-lev-1));\
matrix_type M1_##lev(1<<(lgMtx-lev-1),1<<(lgMtx-lev-1));\
matrix_type M2_##lev(1<<(lgMtx-lev-1),1<<(lgMtx-lev-1));\
matrix_type M3_##lev(1<<(lgMtx-lev-1),1<<(lgMtx-lev-1));\
matrix_type M4_##lev(1<<(lgMtx-lev-1),1<<(lgMtx-lev-1));\
matrix_type M5_##lev(1<<(lgMtx-lev-1),1<<(lgMtx-lev-1));

typedef mat::recursator<matrix_type> recursator_t;

std::vector<recursator_t> M0;
std::vector<recursator_t> M1;
std::vector<recursator_t> M2;
std::vector<recursator_t> M3;
std::vector<recursator_t> M4;
std::vector<recursator_t> M5;
                                                  
allocTempMtx(0);
allocTempMtx(1);
/*allocTempMtx(2);
allocTempMtx(3);
allocTempMtx(4);
allocTempMtx(5);*/

template <typename Matrix>
void print_matrix(Matrix& matrix){ 
  for (int i=0 ; i<matrix.num_rows(); i++ ){
    for(int j=0; j<matrix.num_cols();  j++ ){
      cout.fill (' '); cout.width (8); cout.precision (5); cout.flags (ios_base::left);
      cout << showpoint << /*hex << *reinterpret_cast<long long*>(&*/matrix[i][j] <<"  ";
    }
    cout << endl;
  }
  return;
}


template <typename Matrix>
void fill_matrix_zero(Matrix& matrix){
  for(int i = 0; i < matrix.num_rows(); i++)
    for(int j = 0; j < matrix.num_cols(); j++)
      matrix[i][j] = 0;
}

template <typename Matrix>
void fill_matrix_n(Matrix& matrix){
  for(int i = 0; i < matrix.num_rows(); i++)
    for(int j = 0; j < matrix.num_cols(); j++)
      matrix[i][j] = 2;
}

template < typename Matrix>
void add_base(Matrix & C, Matrix const& A, Matrix const& B)
{
    assert(num_cols(A) == num_cols(B)); assert(num_cols(A) == num_cols(C)); 
    assert(num_rows(A) == num_rows(B)); assert(num_rows(A) == num_rows(C)); 

  for (int i = 0; i < basecasesize && i < num_rows(A); i++)
    for (int j = 0; j < basecasesize && j < num_cols(A); j++)
      C[i][j] = A[i][j] + B[i][j];
}


template <typename Recursator>
void add_m(Recursator C, Recursator A, Recursator B)
{
  if (is_base(A))
  {
    using recursion::base_case_cast;
    typename recursion::base_case_matrix<typename Recursator::matrix_type, BaseTest>::type
      prod= base_case_cast<BaseTest>(C.get_value());
      add_base(prod,
	       base_case_cast<BaseTest>(A.get_value()),
	       base_case_cast<BaseTest>(B.get_value()));
  }
  else 
  {
    add_m(C.north_west(), A.north_west(), B.north_west());
    add_m(C.north_east(), A.north_east(), B.north_east());
    add_m(C.south_east(), A.south_east(), B.south_east());
    add_m(C.south_west(), A.south_west(), B.south_west());
  }
  return;
}

template < typename Matrix>
void sub_base(Matrix & C, Matrix const& A, Matrix const& B)
{
  for (int i = 0; i < basecasesize; i++)
    for (int j = 0; j < basecasesize; j++)
      C[i][j] = A[i][j] - B[i][j];
}

template <typename Recursator>
void sub_m(Recursator C, Recursator A, Recursator B)
{
  if (is_base(A))
  {
    using recursion::base_case_cast;
    typename recursion::base_case_matrix<typename Recursator::matrix_type, BaseTest>::type
      prod= base_case_cast<BaseTest>(C.get_value());
      sub_base(prod,
	       base_case_cast<BaseTest>(A.get_value()),
	       base_case_cast<BaseTest>(B.get_value()));
  }
  else 
  {
    sub_m(C.north_west(), A.north_west(), B.north_west());
    sub_m(C.north_east(), A.north_east(), B.north_east());
    sub_m(C.south_east(), A.south_east(), B.south_east());
    sub_m(C.south_west(), A.south_west(), B.south_west());
  }
  return;  
}



template < typename Matrix>
void strassen_base(Matrix & P, Matrix const& A, Matrix const& B)
{
  fill_matrix_zero(P);
  for (int i = 0; i < basecasesize; i++)
    for (int k = 0; k < basecasesize; k++)     //switching j and k lines cuts
      for (int j = 0; j < basecasesize; j++)   //  run time in half.
	P[i][j] += A[i][k] * B[k][j]; 
}

template <typename Recursator>
void strassen(Recursator P, Recursator A, Recursator B, int lev)
{
  if (A.is_empty() || B.is_empty())
    return;
  if (is_base (A))
  {
    using recursion::base_case_cast;
    typename recursion::base_case_matrix<typename Recursator::matrix_type, BaseTest>::type
      prod= base_case_cast<BaseTest>(P.get_value());
    strassen_base(prod, base_case_cast<BaseTest>(A.get_value()),
		  base_case_cast<BaseTest>(B.get_value()));
    return;
  }
  else
  {
    
    // M1 = (Anw + Ase)*(Bnw + Bse)
    add_m   (M5[lev],        A.north_west(), A.south_east());
    add_m   (M0[lev],        B.north_west(), B.south_east());
    strassen(M1[lev],        M5[lev],        M0[lev],       lev+1);

    // M2 = (Asw + Ase)*Bnw
    add_m   (M5[lev],        A.south_west(), A.south_east());
    strassen(M2[lev],        M5[lev],        B.north_west(), lev+1);

    // M3 = (Asw - Anw)*(Bnw + Bne)
    add_m   (M5[lev],        B.north_west(), B.north_east());
    sub_m   (M0[lev],        A.south_west(), A.north_west());
    strassen(M3[lev],        M0[lev],        M5[lev],       lev+1);

    // M4 = Anw*(Bne - Bse)
    sub_m   (M5[lev],        B.north_east(), B.south_east());
    strassen(M4[lev],        A.north_west(), M5[lev],       lev+1);

    // Pse = M1 - M2 + M4 + M3
    add_m   (M5[lev],        M1[lev],        M4[lev]      );
    add_m   (M0[lev],        M5[lev],        M3[lev]      );
    sub_m   (P.south_east(), M0[lev],        M2[lev]      );

    // M5 = (Anw + Ane)*Bse
    add_m   (M5[lev],        A.north_west(), A.north_east());
    strassen(M3[lev],        M5[lev],        B.south_east(), lev+1);
     
    // Pne = M4 + M5
    add_m   (P.north_east(), M4[lev],        M3[lev]      );

    // M7 = Ase*(Bsw - Bnw)
    sub_m   (M5[lev],        B.south_west(), B.north_west());
    strassen(M4[lev],        A.south_east(), M5[lev],       lev+1);

    // Psw = M2 + M7
    add_m   (P.south_west(), M2[lev],        M4[lev]      );

    // M6 = (Ane - Ase)*(Bsw + Bse)
    add_m   (M5[lev],        B.south_west(), B.south_east());
    sub_m   (M0[lev],        A.north_east(), A.south_east());
    strassen(M2[lev],        M0[lev],        M5[lev],       lev+1);

    // Pnw = M1 + M7 - M5 + M6
    add_m   (M5[lev],        M1[lev],        M4[lev]      );
    add_m   (M0[lev],        M5[lev],        M2[lev]      );
    sub_m   (P.north_west(), M0[lev],        M3[lev]      );
    
    return;
  }
 
}


int test_main(int argc, char* argv[])
{

#define makeRecs(lev)             \
 recursator_t rM0_##lev(M0_##lev);\
 recursator_t rM1_##lev(M1_##lev);\
 recursator_t rM2_##lev(M2_##lev);\
 recursator_t rM3_##lev(M3_##lev);\
 recursator_t rM4_##lev(M4_##lev);\
 recursator_t rM5_##lev(M5_##lev);\
                                  \
 M0.push_back(rM0_##lev);         \
 M1.push_back(rM1_##lev);         \
 M2.push_back(rM2_##lev);         \
 M3.push_back(rM3_##lev);         \
 M4.push_back(rM4_##lev);         \
 M5.push_back(rM5_##lev);          

makeRecs(0);
makeRecs(1);
/*makeRecs(2);
makeRecs(3); 
makeRecs(4);
makeRecs(5);*/


  int order=14; 
  if (argc > 1) order= atoi(argv[1]);

  matrix_type matrix_A(order,order);
  matrix_type matrix_B(order,order); 
  matrix_type matrix_P(order,order);  
  
  fill_matrix_n(matrix_A);
  print_matrix(matrix_A);  printf("\n");
  fill_matrix_n(matrix_B);
  print_matrix(matrix_B);
  fill_matrix_zero(matrix_P);  printf("\n");


  recursion::mat::recursator<matrix_type> recursator_P(matrix_P), 
                                            recursator_A(matrix_A),
                                            recursator_B(matrix_B);

  strassen(recursator_P, recursator_A, recursator_B, 0);

  print_matrix(matrix_P);

  return 0;
}

