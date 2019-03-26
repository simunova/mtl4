// File: simulated_annealing.cpp

#include <iostream>
#include <cmath>
#include <boost/random.hpp>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;
using namespace std;

typedef dense_vector<double>  vector_t;

double inline lg(double x)
{
    return log(x) / log(10.0);
}

double inline square(double x)
{
    return x * x;
}

struct uniform_real
{
    explicit uniform_real(double left, double right)
	: unidist(left, right), uni(rng, unidist)
    {}

    double operator()()
    {
	return uni();
    }

  private:
    boost::mt19937        rng;
    boost::uniform_real<> unidist;
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni;
};



int inline uniform_index()
{
    static boost::mt19937        rng;
    static boost::uniform_int<>  three(0, 2);
    static boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(rng, three);
    return die();
} 


double inline phreeqc_analytic(const vector_t& para, double T)
{
    double lg_k= para[0] + para[1]*T + para[2]/T + para[3]*lg(T) + para[4]/(T*T);
    return pow(10.0, lg_k);
}

double inline min3p_bunsen(const vector_t& para, double T)
{
    double log_beta= para[0] + para[1]*100.0/T + para[2]*log(T/100.0);
    double R= 8.314472 / 101.325; // in L atm / (mol K)
    double Ts= 273.15;
    return exp(log_beta) / (R * Ts);
}

double inline 
error_fun(const vector_t& phreeqc_para, const vector_t& min3p_para, const vector_t& data)
{
    double sum= 0.0;
    for (unsigned i= 0; i < size(data); i++)
	sum+= square(phreeqc_analytic(phreeqc_para, data[i]) - min3p_bunsen(min3p_para, data[i]));
    return sqrt(sum);
}

template <typename Random>
void neighbor_para(vector_t& para, Random random)
{
    para[uniform_index()]+= random();
}

void simulated_annealing(const vector_t& phreeqc_para, vector_t& min3p_para, 
			 const vector_t& data, double dev,
			 double sim_temperature, int max_iter)
{
    double    ref_f= error_fun(phreeqc_para, min3p_para, data), min_f= ref_f;
    vector_t  ref_v(3), min_v(3), new_v(3);
    uniform_real random(-dev, dev);

    int min_counter= 0;
    ref_v= min3p_para;
    for (int i= 0; i < max_iter; i++) {
	new_v= ref_v;
	neighbor_para(new_v, random);
	double new_f= error_fun(phreeqc_para, new_v, data);
	if (new_f < ref_f) {
	    ref_f= new_f; ref_v= new_v;
	    if (new_f < min_f) {
		min_f= new_f; min_v= new_v;
		if (!(min_counter++ % 1000))
		    cout << "New minimum with " << min_f << " at " << new_v;
	    }
	} else 
	    if (rand() < sim_temperature) {
		ref_f= new_f; ref_v= new_v;
	    }
    }
    min3p_para= min_v;
}


int main(int argc, char* argv[])
{
    int max_iter= 1000; 
    if (argc > 1) max_iter= atoi(argv[1]);

    vector_t n2p(5), n2m(3); // N2 parameters for Phreeqc and Min3p Henry coefficients
    n2p[0]= -7.6452, n2p[1]= 7.9606e-3, n2p[2]= n2p[3]= 0, n2p[4]= 1.8604e5;
    n2m[0]= -59.6274, n2m[1]= 85.7661, n2m[2]= 24.3696;  // Values from database
    n2m[0]= -54.0212, n2m[1]= 77.9594, n2m[2]= 21.6595;  // Improved parameters

    vector_t h2p(5), h2m(3); // H2 parameters for Phreeqc and Min3p Henry coefficients
    h2p[0]= -9.3114, h2p[1]= 4.6473e-3, h2p[2]= -4.9335e1, h2p[3]= 1.4341, h2p[4]= 1.2815e+5;
    h2m[0]= -35.3883, h2m[1]= 47.3244, h2m[2]= 14.1778;

    vector_t o2p(5), o2m(3); // O2 parameters for Phreeqc and Min3p Henry coefficients
    o2p[0]= -7.5, o2p[1]= 7.8981e-3, o2p[2]= 0.0, o2p[3]= 0.0, o2p[4]= 2.0075e5;
    o2m[0]= -58.3877, o2m[1]= 85.8079, o2m[2]= 23.8439;

    vector_t temps(201); // Temperatures from 0 to 200 C
    for (unsigned i= 0; i < size(temps); i++)
	temps[i]= 273.15 + double(i);

#if 0
    cout << "Deviation for N2 = " << error_fun(n2p, n2m, temps) << "\n";
    simulated_annealing(n2p, n2m, temps, 0.001, 0.1, max_iter);
    cout << "Deviation for N2 now = " << error_fun(n2p, n2m, temps) 
	 << ", for n2m = " << n2m << "\n";

    cout << "Deviation for H2 = " << error_fun(h2p, h2m, temps) << "\n";
    simulated_annealing(h2p, h2m, temps, 0.001, 0.1, max_iter);
    cout << "Deviation for H2 now = " << error_fun(h2p, h2m, temps) 
	 << ", for h2m = " << h2m << "\n";
#endif

    cout << "Deviation for O2 = " << error_fun(o2p, o2m, temps) << "\n";
    simulated_annealing(o2p, o2m, temps, 0.001, 0.1, max_iter);
    cout << "Deviation for O2 now = " << error_fun(o2p, o2m, temps) 
	 << ", for o2m = " << o2m << "\n";

    return 0;
}

