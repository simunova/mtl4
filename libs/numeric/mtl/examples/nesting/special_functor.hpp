template <>
struct mult_ft<matrix_a_type, matrix_b_type, matrix_c_type>
{
    void operator()(matrix_a_type const& a, matrix_b_type const& b, matrix_c_type& c) 
    { /* Faster code for this type triplet ... */ }
};
