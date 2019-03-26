#include <iostream>
#include <fstream>
#include <utility>
#include <sstream>

#include <boost/numeric/mtl/mtl.hpp>

typedef double value_type;
typedef std::size_t size_type;
typedef mtl::dense_vector<value_type>  vector_type;
typedef mtl::compressed2D<value_type>  matrix_type;
const double pi(3.14159265358979323846);

template <typename Vector, typename Matrix>
class grad_f_ftor
{
  public:
    grad_f_ftor(const Matrix& M, const Matrix& K, const Matrix& G, const value_type& beta, const value_type& w) 
      : M(M), K(K), G(G), beta(beta), w(w)
      {
	q.change_dim(num_cols(G));
	x.change_dim(num_cols(G));
	x=0.0; q=0.0;
	for(size_type i=1; i< size(x); i++)
	  x(i)= x(i-1) + 0.01;
      }
      
    //get timedependent input q for the heating
    Vector get_input(const value_type& time) 
    {
      value_type center(0.5*sin(2*pi*time)+0.5);
      for(size_type i=0; i< size(x); i++)
	  q(i)=beta*cos(std::min(std::abs(pi*(x(i)-center)/(2*w)),pi/2));
      return q;
    }
    template <typename VectorIn>
    Vector operator()(const VectorIn& u, const value_type& time) 
    {
	q=get_input(time);
	Vector x(K*u + G*q);
	return x;
    }
 
  private:
    Vector      q, x;
    Matrix      M, K, G;
    value_type  beta, w;
   
};


template <typename grad_f_ftor, typename Vector, typename Matrix>
Vector ode23s(grad_f_ftor func, value_type start_time, value_type end_time, Vector start_value, Matrix M, Matrix K){

  value_type time(start_time), time_step(0.005), h= 0.005, gamma(1-1/sqrt(2));
  mtl::dense2D<value_type> LU(M-h*gamma*K);
  lu(LU);
  Vector x(start_value);
  size_type time_counter(0);
  Vector k1, k2, step1, step2;
  while (time < end_time){
      k1= lu_solve_apply(LU, func(x,time));
      step1= x + h * k1;
      step2= func(step1, time + time_step) - 2*M*k1;
      k2= lu_solve_apply(LU, step2);
      save_data(x, time_counter);
      x+= 3/2*h*k1 + 1/2*h*k2; 
      time_counter++;
      time+= time_step;
  }
  return x;
}

//save current state of solution x(time)
template<typename Vector>
void save_data(const Vector& x, const size_type& time) 
{
    std::stringstream name;
    name  <<"plot_data_" << time <<".dat";
#if 1
    std::cout<< "#!/usr/bin/gnuplot\n";
    std::cout<< "set term png\n";
   // std::cout<< "set zrange[-4:14]\n";
    std::cout<< "set output \""<<time<<".png\"\n";
    std::cout<< "splot './plot_data_"<< time<<".dat' using 1:2:3 with pm3d\n";
#endif
    std::fstream f;
    f.open(name.str().c_str(), std::ios::out);
    size_type n(1111), n1(101);
 
    size_type row(0), col(0);
    for(size_type i=0; i < n; i++){
      if(i%(n/n1)==0){
	row=0;col+=1;
	f << "\n";
      }
      f <<  col << " " << row << " " << x(i) << "\n";
      row+=1;
    }
    f.close();   
}

 
int main( int  , char ** )
{
     //read jörgs matrices
    matrix_type M(mtl::io::matrix_market( "M.mtx")),
		K(mtl::io::matrix_market( "K.mtx")),
		G(mtl::io::matrix_market( "G.mtx"));
   
    size_type n(num_rows(K));
    value_type  start_time(0.0), end_time(1.0);
    vector_type  x(n,0.0), x0(n,0.0);
  
    grad_f_ftor< vector_type, matrix_type >  grad_f(M, K, G, 0.2, 0.1);
    
    x= ode23s(grad_f, start_time, end_time, x0, M, K);
    
    return 0;
}