#include <iostream>
#include <cmath>

using namespace std;

template <typename T>

class polynom{

   public:
    T* data;
    explicit polynom(int size) : size(size), data(new T[100]) {}

    polynom(const polynom& p) : size(p.size), data(new T[100]) 
    {
	for(int i= 0; i <= size; i++)
	    data[i]= p.data[i];
    }

    double get_size() const { return size; }

    double set_size(int s){ size=s; }

    polynom stamm(){
	polynom b(size+1);
	for(int i=1; i<=b.get_size();i++){
		b.data[i]=data[i-1]/i;
	}
	return b;
    }

    T integral(T xu, T xo){
	T y=T(0.0);
	polynom b=stamm();
	for(int i=0;i<=b.get_size();i++){
		y+=b.data[i]*(pow(xo,i)-pow(xu,i));
	}
	return y;
    }
    
    polynom diff(){
	polynom b(size-1);
	for(int i=0; i<=b.get_size();i++){
		b.data[i]=(i+1)*data[i+1];
	}
	return b;
    }

    T wert(T a){
	T y=T(0.0);
	for(int i=0;i<size;i++){
		y+=data[i]*(pow(a,i));
	}
	return y;
    }

    polynom operator+(const polynom& p) const
    {
	int max;
	if(size>p.get_size()){
		max=size;
	}else{
		max=p.get_size();
	}
	polynom c(max);
	for(int i=0; i<=max; i++){
		c.data[i]=data[i]+p.data[i];
	}
	return c;
	
    }

    polynom operator*(const polynom& p) const
    {
	polynom c(size+p.get_size());
	for (int i=0; i<= size+p.get_size(); i++){
		for( int j=0; j<=i; j++){
			c.data[i]+=(data[j])*(p.data[i-j]);
		}
	}
	return c;
    }

    polynom operator*(T b) const
    {
	polynom c(size);
	for(int i=0;i<=c.get_size();i++){
		c.data[i]=data[i]*b;
	}
	return c;
    }


  private:
    int size;
    
};


