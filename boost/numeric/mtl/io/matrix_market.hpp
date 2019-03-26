// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#ifndef MTL_IO_MATRIX_MARKET_INCLUDE
#define MTL_IO_MATRIX_MARKET_INCLUDE

#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <locale>
#include <complex>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
// #include <boost/algorithm/string/case_conv.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <boost/numeric/mtl/io/matrix_file.hpp>
#include <boost/numeric/mtl/io/read_filter.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/string_to_enum.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>

namespace mtl { namespace io {


/// Input file stream for files in matrix market format
class matrix_market_istream
{
    class pattern_type {};
    typedef matrix_market_istream        self;
    void check_stream(const std::string& file_name= std::string()) // not const to delete new_stream
    {
	if (!my_stream.good()) {
	    std::string message("matrix_market_istream: Error in input stream!\n");
	    if (file_name != std::string())
		message+= "Probably file " + file_name + " not found.\n";
	    std::cerr << message;
	    if (new_stream)
		delete new_stream, new_stream= 0; // To get valgrind quite
	    throw(file_not_found(message.c_str()));
	}
    }

  public:
    explicit matrix_market_istream(const char* p) : new_stream(new std::ifstream(p)), my_stream(*new_stream) { check_stream(p); }
    explicit matrix_market_istream(const std::string& s) : new_stream(new std::ifstream(s.c_str())), my_stream(*new_stream) { check_stream(s); }
    explicit matrix_market_istream(std::istream& s= std::cin) : new_stream(0), my_stream(s) { check_stream(); }

    ~matrix_market_istream() 
    { 	if (new_stream) delete new_stream;    }

    template <typename Coll>
    self& operator>>(Coll& c) 
    { 	return read(c, typename mtl::traits::category<Coll>::type());    }

    /// Close only my own file, i.e. if filename and not stream is passed in constructor
    void close() { if (new_stream) new_stream->close(); }

  protected:
    template <typename Matrix> self& read(Matrix& A, tag::matrix);

    void to_lower(std::string& s) const
    {
	using namespace std;
	for (unsigned i= 0; i < s.size(); i++)
	    s[i]= tolower(s[i]);
    }

    void set_symmetry(std::string& symmetry_text)
    {
	to_lower(symmetry_text); 
	const char* symmetry_options[]= {"general", "symmetric", "skew-symmetric", "hermitian"};
	my_symmetry= string_to_enum(symmetry_text, symmetry_options, symmetry());
    }

    void set_sparsity(std::string& sparsity_text)
    {
	to_lower(sparsity_text); 
	const char* sparsity_options[]= {"coordinate", "array"};
	my_sparsity= string_to_enum(sparsity_text, sparsity_options, sparsity());
    }

    template <typename Inserter, typename Value>
    void read_matrix(Inserter& ins, Value)
    {
	typedef typename Collection<typename Inserter::matrix_type>::size_type size_type;
	read_filter<Inserter> filter(ins);
	if (my_sparsity == coordinate) // sparse
	    // while (my_stream) { // sometimes does an extra erroneous loop
	    for (std::size_t i= 0; i < nnz; i++) {
		size_type r, c;
		my_stream >> r >> c;
		if (!my_stream) break; // in case while(my_stream) caught an empty line at the end
		insert_value(ins, r-1, c-1, filter, Value());
	    }
	else // dense 
	    for (std::size_t c= 0; c < ncols; c++) 
		for (std::size_t r= 0; r < nrows; r++)
		    insert_value(ins, r, c, filter, Value());
    }

    template <typename Inserter, typename Filter, typename Value>
    void insert_value(Inserter& ins, std::size_t r, std::size_t c, const Filter& filter, Value) 
    {
	using mtl::conj;
	typedef typename Collection<typename Inserter::matrix_type>::value_type mvt;
	Value v;
	read_value(v);
	// std::cout << "Going to insert at [" << r << "][" << c << "] value " << which_value(v, mvt()) << "\n";
	if (filter(r, c))
	    ins[r][c] << which_value(v, mvt());
	if (r != c && filter(c, r)) 
	    switch(my_symmetry) {
	      case symmetric:      ins[c][r] << which_value(v, mvt()); break;
	      case skew:           ins[c][r] << -which_value(v, mvt()); break;
	      case Hermitian:      ins[c][r] << conj(which_value(v, mvt())); break;
	      default:             ; // do nothing
	    }
    }

    void read_value(pattern_type) {}
    void read_value(double& v) { my_stream >> v;}
    void read_value(long& v) { my_stream >> v;}
    void read_value(std::complex<double>& v) 
    { 
	double r, i; my_stream >> r >> i; v= std::complex<double>(r, i);
    }

    // Which value to be inserted? Itself if exist and 0 for pattern; complex are 
    template <typename Value, typename MValue> MValue which_value(Value v, MValue) { return boost::numeric_cast<MValue>(v); }
    template <typename MValue> MValue which_value(pattern_type, MValue) { return boost::numeric_cast<MValue>(0.0); }
    template <typename MValue> MValue which_value(std::complex<double>, MValue) { MTL_THROW(runtime_error("Cannot convert complex value in real\n")); return 1; }
    std::complex<long double> which_value(std::complex<double> v, std::complex<long double>) { return boost::numeric_cast<std::complex<long double> >(v); }
    std::complex<double> which_value(std::complex<double> v, std::complex<double>) { return v; }
    std::complex<float> which_value(std::complex<double> v, std::complex<float>) { return std::complex<float>(float(real(v)), float(imag(v))); }

    std::ifstream      *new_stream;
    std::istream       &my_stream;
    enum symmetry {general, symmetric, skew, Hermitian} my_symmetry;
    enum sparsity {coordinate, array} my_sparsity;
    std::size_t nrows, ncols, nnz;
};



// Matrix version
template <typename Matrix>
matrix_market_istream& matrix_market_istream::read(Matrix& A, tag::matrix)
{
    std::string marker, type, sparsity_text, value_format, symmetry_text;
    my_stream >> marker >> type >> sparsity_text >> value_format >> symmetry_text;
#if 0    
    std::cout << marker << ", " << type << ", " << sparsity_text << ", " 
	      << value_format << ", " << symmetry_text << "\n";
#endif
    MTL_THROW_IF(marker != std::string("%%MatrixMarket"), 
		 runtime_error("File not in Matrix Market format"));
    MTL_THROW_IF(type != std::string("matrix"), 
		 runtime_error("Try to read matrix from non-matrix file"));

    set_symmetry(symmetry_text);
    set_sparsity(sparsity_text);

    char first, comment[80];
    do {
	my_stream >> first;
	if (first == '%') // comments start with % -> ignore them
	    my_stream.getline(comment, 80, '\n'); // read rest of line
	else
	    my_stream.putback(first); // if not commment we still need it
    } while (first == '%');

    my_stream >> nrows >> ncols;
    // std::cout << nrows << "x" << ncols << ", " << nnz << " non-zeros\n";	
    A.change_dim(nrows, ncols);
    set_to_zero(A);

    std::size_t slot_size;
    if (sparsity_text == std::string("coordinate")) {
	my_stream >> nnz; slot_size= std::max(std::size_t(double(nnz) / double(A.dim1()) * 1.25), std::size_t(1));
    } else
	slot_size= A.dim2(); // maximal value (if A is dense it does not matter anyway)

    // Create enough space in sparse matrices
    mat::inserter<Matrix> ins(A, slot_size);

    if (value_format == std::string("real"))
	read_matrix(ins, double());
    else if (value_format == std::string("integer"))
	read_matrix(ins, long());
    else if (value_format == std::string("complex"))
	read_matrix(ins, std::complex<double>());
    else if (value_format == std::string("pattern"))
	read_matrix(ins, pattern_type());
    else
	MTL_THROW(runtime_error("Unknown tag for matrix value type in file"));

    return *this;
}


class matrix_market_ostream 
{
    typedef matrix_market_ostream        self;
public:
    explicit matrix_market_ostream(const char* p) : new_stream(new std::ofstream(p)), my_stream(*new_stream) {}
    explicit matrix_market_ostream(const std::string& s) : new_stream(new std::ofstream(s.c_str())), my_stream(*new_stream) {}
    explicit matrix_market_ostream(std::ostream& s= std::cout) : new_stream(0), my_stream(s) {}

    ~matrix_market_ostream() { if (new_stream) delete new_stream; }

    template <typename Coll>
    self& operator<<(const Coll& c) 
    { 
	return write(c, typename mtl::traits::category<Coll>::type());
    }

    /// Close only my own file, i.e. if filename and not stream is passed in constructor
    void close() { if (new_stream) new_stream->close(); }

private:
    template <typename Matrix> self& write(const Matrix& A, tag::matrix)
    {
	matrix_status_line(A);
	if (sparsity(A) == std::string("coordinate "))
	    return write_sparse_matrix(A);
	else
	    return write_dense_matrix(A);
    }

    template <typename Matrix> self& write_sparse_matrix(const Matrix& A)
    {
	my_stream << num_rows(A) << " " << num_cols(A) << " " << A.nnz() << "\n";
	
	typename mtl::traits::row<Matrix>::type             row(A); 
	typename mtl::traits::col<Matrix>::type             col(A); 
	typename mtl::traits::const_value<Matrix>::type     value(A); 
	typedef typename mtl::traits::range_generator<tag::major, Matrix>::type  cursor_type;
	typedef typename mtl::traits::range_generator<tag::nz, cursor_type>::type icursor_type;

	for (cursor_type cursor = begin<tag::major>(A), cend = end<tag::major>(A); cursor != cend; ++cursor)
	    for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor)
		my_stream << row(*icursor)+1 << " " << col(*icursor)+1 << " ", write_value(value(*icursor)), my_stream << "\n";
	return *this;
    }

    template <typename Matrix> self& write_dense_matrix(const Matrix& A)
    {
	my_stream << num_rows(A) << " " << num_cols(A) << "\n";
	for (std::size_t c = 0; c < num_cols(A); ++c)
	    for (std::size_t r = 0; r < num_rows(A); ++r) {
		write_value(A[r][c]), my_stream << " ";
	    my_stream << "\n";
	}
	return *this;
    }

    template <typename Value>
    typename boost::enable_if<boost::is_integral<Value> >::type
    write_value(const Value& v) { my_stream << v; }

    template <typename Value>
    typename boost::enable_if<boost::is_floating_point<Value> >::type
    write_value(const Value& v) 
    { 
	my_stream.precision(std::numeric_limits<Value>::digits10 + 1);
	my_stream.setf(std::ios::scientific); 
	my_stream << v; 
	my_stream.unsetf(std::ios::scientific); 
    }

    template <typename Value>
    void write_value(const std::complex<Value>& v) 
    { 
	my_stream.precision(std::numeric_limits<Value>::digits10 + 1);
	my_stream.setf(std::ios::scientific); 
	my_stream << real(v) << " " << imag(v); 
	my_stream.unsetf(std::ios::scientific); 
    }

    // Will be generalized via traits::is_symmetric and alike
    template <typename Matrix>
    std::string symmetry(const Matrix&) const { return std::string("general\n"); }

    template <typename Matrix>
    std::string sparsity(const Matrix&) const 
    {
	return std::string( mtl::traits::is_sparse<Matrix>::value ? "coordinate " : "array " );
    }

    template <typename Value>
    typename boost::enable_if<boost::is_integral<Value>, std::string>::type
    value(const Value&) const { return std::string("integer "); }

    template <typename Value>
    typename boost::enable_if<boost::is_floating_point<Value>, std::string>::type
    value(const Value&) const { return std::string("real "); }

    template <typename Value>
    std::string value(const std::complex<Value>&) const { return std::string("complex "); }

    template <typename Matrix>
    void matrix_status_line(const Matrix& A) const
    {
	typedef typename Collection<Matrix>::value_type value_type;
	std::string st(std::string("%%MatrixMarket  matrix ") + sparsity(A) + value(value_type()) + symmetry(A));
	my_stream << st;
    }

protected:
    std::ofstream      *new_stream;
    std::ostream       &my_stream;
};


}} // namespace mtl::io

#endif // MTL_IO_MATRIX_MARKET_INCLUDE
