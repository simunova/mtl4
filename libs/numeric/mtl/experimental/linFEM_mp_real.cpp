#include <iostream>
#include <cmath>
#include <cassert>
#include <arprec/mp_real.h>
#include "Polynomklasse_templ.hpp"
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;

class mp_starter
{
  public:
    mp_starter() { mp::mp_init(40); }
    ~mp_starter() { mp::mp_finalize(); }
};

const mp_starter dummy;

int const gp=10; //=N+1
mp_real gitter[gp];
mp_real e("1.0");


polynom<mp_real> phi(int i, mp_real eins, mp_real zwei){
	polynom<mp_real> c(1);
	if(i==0||i==gp-1){
		c.data[0]=0.0;
		c.data[1]=0.0;
		return c;
	}
	if(eins==gitter[i] && zwei==gitter[i+1]){
		c.data[0]=-gitter[i+1]/(gitter[i]-gitter[i+1]);
		c.data[1]=1/(gitter[i]-gitter[i+1]);
	}else if(eins==gitter[i-1] && zwei==gitter[i]){
		c.data[0]=-gitter[i-1]/(gitter[i]-gitter[i-1]);
		c.data[1]=1/(gitter[i]-gitter[i-1]);
	}else{
		c.data[0]=0.0;
		c.data[1]=0.0;
	}
	return c;
}

mp_real func_c(mp_real x){
	return 1.0;
}

mp_real f(mp_real x){
	return 1.0;
}

polynom<mp_real> inter_c(int i, mp_real eins, mp_real zwei){
	polynom<mp_real> c(1);
	if(i==0||i==gp-1){
		c.data[0]=0.0;
		c.data[1]=0.0;
		return c;
	} 
	if(eins==gitter[i] && zwei==gitter[i+1]){
		c.data[0]=func_c(gitter[i])*(-1.0)*(gitter[i+1])/(gitter[i]-gitter[i+1]);
		c.data[1]=func_c(gitter[i])/(gitter[i]-gitter[i+1]);
	}else if(eins==gitter[i-1] && zwei==gitter[i]){
		c.data[0]=func_c(gitter[i])*(-1.0)*(gitter[i-1])/(gitter[i]-gitter[i-1]);
		c.data[1]=func_c(gitter[i])/(gitter[i]-gitter[i-1]);
	}else{
		c.data[0]=0.0;
		c.data[1]=0.0;
	}
	return c;
}
polynom<mp_real> inter_f(mp_real eins, mp_real zwei){
	polynom<mp_real> c(1);
	c.data[0]=f(eins)-(f(zwei)-f(eins))/(zwei-eins)*eins;
	c.data[1]=(f(zwei)-f(eins))/(zwei-eins);
	return c;

}



mp_real eintrag(int i, int j){
	mp_real tmp=0.0;
	for(int k=0;k<gp-1;k++){
		polynom<mp_real> da=phi(i,gitter[k],gitter[k+1]).diff();
		polynom<mp_real> db=phi(j,gitter[k],gitter[k+1]).diff();
		polynom<mp_real> b=phi(j,gitter[k],gitter[k+1]);
		polynom<mp_real> I_ca=inter_c(i,gitter[k],gitter[k+1]);
		tmp+=(e*e)*(da*db).integral(gitter[k],gitter[k+1])+(I_ca*b).integral(gitter[k],gitter[k+1]);
	}
	return tmp;	
	
}

mp_real vektoreintrag(int j){
	mp_real tmp=0.0;
	for(int k=0;k<gp-1;k++){
		polynom<mp_real> b=phi(j,gitter[k],gitter[k+1]);
		polynom<mp_real> I_f=inter_f(gitter[k],gitter[k+1]);
		tmp+=(I_f*b).integral(gitter[k],gitter[k+1]);
	}
	return tmp;
}

int main(){
	cout<<"test";
	// mp::mp_init(40);
	cout<<"test1";
	for (int i=0;i<gp;i++){
		cout<<"test2";
		gitter[i]=1.0*i/(gp-1);
		cout<<"test3";
		cout<<"gitter["<<i<<"]="<<gitter[i]<<" \n";
	
	}

	cout.precision(3);	
	
	
	mtl::dense2D<mp_real> matrix(gp-2,gp-2);

	for(int i=1;i<gp-1;i++){
		for(int j=1;j<gp-1;j++){
			matrix[i-1][j-1]=eintrag(i,j);
			
		}
	}
	
	cout<<matrix;
	
	mtl::dense_vector<mp_real> vector(gp-2),p(gp-2);
	
	for(int i=1;i<gp-1;i++){
		vector[i-1]=vektoreintrag(i);
	}
	cout<<vector<<"\n";
	
	p=lu_solve(matrix,vector);
	
	cout<<p<<"\n";

	
	return 0;

	
}
