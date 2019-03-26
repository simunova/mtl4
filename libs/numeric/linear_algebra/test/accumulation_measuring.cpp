
#include <libs/numeric/linear_algebra/test/accumulation.hpp>

#include <boost/timer.hpp>
#include <numeric>
#include <complex>
#include <fstream>


// The following code is for measuring
// ######################################################################

template <typename Element>
struct magma_add : public math::add<Element> {};

template <typename Element>
struct magma_mult : public math::mult<Element> {};

template <typename Element>
struct magma_max : public math::max<Element> {};


// I really shouldn't need this
namespace math {
    template <typename Element>
    struct identity_t<magma_add<Element>, Element> 
	: public identity_t<add<Element>, Element> {};

    template <typename Element>
    struct identity_t<magma_mult<Element>, Element> 
	: public identity_t<mult<Element>, Element> {};

    template <typename Element>
    struct identity_t<magma_max<Element>, Element> 
	: public identity_t<max<Element>, Element> {};
}

template <typename T> struct magma_type {};
template <typename T> struct magma_type<math::add<T> > {  typedef magma_add<T> type;};
template <typename T> struct magma_type<math::mult<T> > {  typedef magma_mult<T> type;};
template <typename T> struct magma_type<math::max<T> > {  typedef magma_max<T> type;};

namespace mtl {


// Dispatching between simple and unrolled version
template <typename Iter, typename Value, typename Op>
  where std::ForwardIterator<Iter> 
                  && std::Convertible<Value, std::ForwardIterator<Iter>::value_type>
                  && math::Magma<Op, std::ForwardIterator<Iter>::value_type>
typename std::ForwardIterator<Iter>::value_type 
inline fast_accumulate(Iter first, Iter last, Value init, Op op)
{
    // std::cout << "Simple accumulate\n";
    return mtl::accumulate_simple(first, last, init, op);
}

  // The last 2 are for resolving ambiguities
template <typename Iter, typename Value, typename Op>
    where  std::RandomAccessIterator<Iter> 
	          && std::Convertible<Value, std::RandomAccessIterator<Iter>::value_type>
		  && math::CommutativeMonoid<Op, std::RandomAccessIterator<Iter>::value_type> 
typename std::RandomAccessIterator<Iter>::value_type 
inline fast_accumulate(Iter first, Iter last, Value init, Op op)
{
    // std::cout << "Unrolled accumulate\n";
    return mtl::accumulate_unrolled(first, last, init, op);
}

} // namespace mtl

template<typename OStream, typename Operation, typename Element> 
Element timing_accumulation(OStream& os, Element init, Operation op, int vector_size, int repetitions)
{
    std::vector<Element>   v(vector_size);
    boost::timer start;
    Element result;
    for (int i= 0; i < repetitions; i++)
	result= mtl::fast_accumulate(v.begin(), v.end(), init, Operation());
    double duration = start.elapsed();
    os << duration / repetitions * 1000000; // << "µs" << "\n";
    return result;
}

template<typename OStream, typename Operation, typename Element> 
Element timing_stl_accumulation(OStream& os, Element init, Operation op, int vector_size, int repetitions)
{
    std::vector<Element>   v(vector_size);
    boost::timer start;
    Element result;
    for (int i= 0; i < repetitions; i++)
	result= std::accumulate(v.begin(), v.end(), init, Operation());
    double duration = start.elapsed();
    os << duration / repetitions * 1000000; // << "µs" << "\n";
    return result;
}


template<typename OStream, typename Operation, typename Element> 
Element timing_all_dispatchings(OStream& os, Element init, Operation op, int vector_size, int repetitions)
{
    typedef typename magma_type<Operation>::type   magma_op;

    os << vector_size << ", ";
    timing_accumulation(os, init, op, vector_size, repetitions);  os << ", ";
    timing_accumulation(os, init, magma_op(), vector_size, repetitions);  os << ", ";
    timing_stl_accumulation(os, init, magma_op(), vector_size, repetitions);  os << "\n";
}


struct ostream_repeater
{
    std::ostream& os1, &os2;
    ostream_repeater(std::ostream& os1, std::ostream& os2) : os1(os1), os2(os2) {}
};

template <typename T>
ostream_repeater& operator<< (ostream_repeater& d, T v)
{
    d.os1 << v; d.os2 << v; return d;
}


template<typename Operation, typename Element> 
Element timing_series(Element init, Operation op, int max_size, int step_size, int repetitions,
		      const char* file_name)
{
    std::cout << "Will write in file: " << file_name << "\n";
    std::ofstream myfile;
    myfile.open (file_name);

    ostream_repeater  repeater(std::cout, myfile);

    for (int i= step_size; i <= max_size; i+= step_size)
	timing_all_dispatchings(repeater, init, op, i, repetitions);
    myfile.close();

}



int main(int argc, char* argv[])
{
    using std::cout;
    using math::identity;

    if (argc < 4) {cout << "usage: accumulation_measuring <# of measurements> <max_size> <step_size>\n"; exit(1); }
    int repetitions= atoi(argv[1]), max_size= atoi(argv[2]), step_size= atoi(argv[3]);

    typedef std::complex<double> c_t;
    std::complex<double>  c0, c1(1.0);
    // timing_all_dispatchings(cout, c0, math::add<std::complex<double> >(), 1000, repetitions);
    // timing_all_dispatchings(cout, c1, math::mult<std::complex<double> >(), 1000, repetitions);

    timing_series(0.0f, math::add<float>(), max_size, step_size, repetitions, "add_float.dat");
    timing_series(1.0f, math::mult<float>(), max_size, step_size, repetitions, "mult_float.dat");
    typedef math::max<float>  mf_t;
    timing_series(identity(mf_t(), 0.0f), mf_t(), max_size, step_size, repetitions, "max_float.dat");

    timing_series(0.0, math::add<double>(), max_size, step_size, repetitions, "add_double.dat");
    timing_series(1.0, math::mult<double>(), max_size, step_size, repetitions, "mult_double.dat");
    //timing_series(0.0f, math::max<double>(), max_size, step_size, repetitions, "max_double.dat");

    timing_series(0, math::add<int>(), max_size, step_size, repetitions, "add_int.dat");
    timing_series(1, math::mult<int>(), max_size, step_size, repetitions, "mult_int.dat");
    typedef math::max<int>  mi_t;
    // timing_series(identity(mi_t(), int(0)), mi_t(), max_size, step_size, repetitions, "max_int.dat");

    timing_series(c0, math::add<c_t>(), max_size, step_size, repetitions, "add_complex_double.dat");
    timing_series(c1, math::mult<c_t>(), max_size, step_size, repetitions, "mult_complex_double.dat");


    return 0;
}
