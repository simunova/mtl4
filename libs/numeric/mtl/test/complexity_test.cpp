#include <iostream>
#include <stdio.h>

#include <boost/numeric/mtl/utility/complexity.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/less.hpp>


using namespace std;
using namespace mtl::complexity_classes;
namespace mpl = boost::mpl;

struct wahr {
  void operator() () {
    cout << "vrai\n"; }
};

struct falsch {
  void operator() () {
    cout << "faute\n"; }
};

void schreib(infinite) {
  cout << "unendlich\n"; }

void schreib(polynomial) {
    cout << "polynomial\n"; }

void schreib(quadratic) {
    cout << "quadratisch\n"; }

void schreib(n_polylog_n) {
    cout << "n log^k n\n"; }
 
void schreib(n_log_n) {
    cout << "n log n\n"; }
 
void schreib(linear) {
    cout << "linear\n"; } 

void schreib(linear_cached) {
    cout << "linear cached\n"; }

void schreib(polylog_n) {
    cout << "polynomial log\n"; }

void schreib(log_n) {
    cout << "log\n"; }

void schreib(constant) {
    cout << "constant\n"; }

void schreib(cached) {
    cout << "cached\n"; }



template <typename X, typename Y> void write_less(X, Y) {
  typedef mpl::less<X, Y> less_res;
  typename mpl::if_<less_res, wahr, falsch>::type()();
}

template <typename X, typename Y> void write_plus(X, Y) {
  schreib(typename mtl::complexity_classes::plus<X, Y>::type());
}

template <typename X, typename Y> void write_mal(X, Y) {
  schreib(typename mtl::complexity_classes::times<X, Y>::type());
}

 
int main (int, char**) {

  write_less(quadratic(), infinite());
  write_less(quadratic(), linear());

  write_plus(quadratic(), infinite());
  write_plus(quadratic(), linear());
  write_plus(linear(), quadratic());
  write_plus(n_log_n(), linear());
  
  write_mal(quadratic(), infinite());
  write_mal(infinite(), quadratic());
  write_mal(quadratic(), quadratic());
  write_mal(linear(), quadratic());
  write_mal(linear(), log_n());
  write_mal(linear(), polylog_n());
  write_mal(n_log_n(), log_n());
  write_mal(n_log_n(), polylog_n());
  write_mal(cached(), log_n());
  write_mal(log_n(), log_n());
  write_mal(cached(), constant());



  return 0;
}
