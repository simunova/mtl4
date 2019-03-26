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
//
// This example is written by Andreas Byfut.

# include <array>
# include <vector>

# include <boost/numeric/itl/itl.hpp>
# include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;
using namespace itl;

typedef long t_int;
typedef long double t_double;

typedef dense_vector<t_double,parameters<tag::col_major, fixed::dimension<2>>> t_vector_2;
typedef dense2D<t_double,mat::parameters<tag::row_major,mtl::index::c_index,mtl::fixed::dimensions<3,3>>> t_matrix_3x3;
typedef dense2D<t_double,mat::parameters<tag::row_major,mtl::index::c_index,mtl::fixed::dimensions<3,2>>> t_matrix_3x2;

typedef std::vector<t_vector_2> t_c4n;
typedef std::vector<std::array<t_int,3>> t_n4e;
typedef std::vector<t_double> t_u4n;
typedef std::vector<std::array<t_int,2>> t_n4ed;


void write_svg_color( std::array<t_int,3>& colorVal, const t_double& percentage )
{
    if ( percentage < 0.5 )
    {
        colorVal[0] = 2*percentage * 255;
        colorVal[1] = 2*percentage * 255;
        colorVal[2] = 255;
    }
    else
    {
        colorVal[0] = 255;
        colorVal[1] = 2*(1-percentage) * 255;
        colorVal[2] = 2*(1-percentage) * 255;
    }
}

void write_svg( const t_c4n& c4n, const t_n4e& n4e, const t_u4n& u4n, const std::string& filename )
{
    t_double scaling = 500;
    t_double minX = c4n[0](0), maxX = c4n[0](0), minY = c4n[0](1), maxY = c4n[0](1);
    t_double minU = u4n[0], maxU = u4n[0];
    for ( t_int curNode = 1, nrNodes = c4n.size(); curNode < nrNodes; ++curNode )
    {
        if ( c4n[curNode](0) < minX ) minX = c4n[curNode](0);
        if ( c4n[curNode](0) > maxX ) maxX = c4n[curNode](0);
        if ( c4n[curNode](1) < minY ) minY = c4n[curNode](1);
        if ( c4n[curNode](1) > maxY ) maxY = c4n[curNode](1);
        if ( u4n[curNode] < minU ) minU = u4n[curNode];
        if ( u4n[curNode] > maxU ) maxU = u4n[curNode];
    }
    std::cout << "minU " << minU << "\n";
    std::cout << "maxU " << maxU << "\n";

    std::ofstream xmlFile(filename.c_str());
    if ( xmlFile.is_open() )
    {
        xmlFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
        xmlFile << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.0\" width=\"" << scaling*(maxX-minX)+100 << "\" height=\"" << scaling*(maxY-minY) << "\">\n";

        // Draw individual elements colored with the solution at the midpoint.
        for ( t_int curElem = 0, nrElem = n4e.size(); curElem < nrElem; ++curElem )
        {
            t_double u4e = (u4n[n4e[curElem][0]]+u4n[n4e[curElem][1]]+u4n[n4e[curElem][2]])/3.0;
            t_double u4e_percent = (u4e-minU)/(maxU-minU);
            std::array<t_int,3> colorVal; write_svg_color(colorVal,u4e_percent);
            xmlFile << "<polygon points=\"";
            xmlFile << scaling*(c4n[n4e[curElem][0]](0)-minX) << ", " << scaling*(maxY-c4n[n4e[curElem][0]](1)) << " ";
            xmlFile << scaling*(c4n[n4e[curElem][1]](0)-minX) << ", " << scaling*(maxY-c4n[n4e[curElem][1]](1)) << " ";
            xmlFile << scaling*(c4n[n4e[curElem][2]](0)-minX) << ", " << scaling*(maxY-c4n[n4e[curElem][2]](1));
            xmlFile << "\" style=\"fill:rgb(" << colorVal[0] << "," << colorVal[1] << "," << colorVal[2] << ");stroke:black;stroke-width:1;stroke-linecap='round'\"/>\n";
        }

        // Draw legend for the solution.
        for ( t_int i = 0, I = 15; i <= I; ++i )
        {
            std::array<t_int,3> colorVal; write_svg_color(colorVal,t_double(i)/t_double(I));
            xmlFile << "<rect x=\"" << scaling*(maxX-minX)+20 << "\" y=\"" << scaling*(maxY-minY-(i+1)*(maxY-minY)/(I+1)) << "\" width=\"80\" height=\"" << scaling*(maxY-minY)/(I+1) <<
                       "\" style=\"fill:rgb(" << colorVal[0] << "," << colorVal[1] << "," << colorVal[2] << ");stroke:black;stroke-width:1;stroke-linecap='round'\"/>";
            xmlFile << "<text x=\"" << scaling*(maxX-minX)+60 << "\" y=\"" << scaling*(maxY-minY-(i+0.45)*(maxY-minY)/(I+1)) << "\" text-anchor=\"middle\" font-size=\"12\">" << (maxU-minU)*i/I+minU << "</text>\n";
        }

        xmlFile << "</svg>";
        xmlFile.close();
    }
}

//! Make a red-refinment of the given mesh
/*! The refinement is depicted below:
     /\                      /\
    /  \                    /  \
   /    \    ------->      /----\
  /      \                / \  / \
 /        \              /   \/   \
------------            ------------ */
void red_refine( t_c4n& c4n, t_n4e& n4e, t_n4ed& dbnd, t_n4ed& nbnd )
{
    // Determine all edges of the mesh.
    compressed2D<t_int> K(c4n.size(),c4n.size()), T;
    {
        mat::inserter<compressed2D<t_int>> ins(K);
        for ( t_int curElem = 0, nrElem = n4e.size(); curElem < nrElem; ++curElem )
            for ( t_int i = 0; i < 3; ++i )
                ins[n4e[curElem][i]][n4e[curElem][(i+1)%3]] << 1;
    }
    T = K + trans(K);
    K = upper(T);
    // Compute mid-points for each edge and save indices.
    c4n.reserve( c4n.size() + K.nnz() );
    typedef typename mtl::traits::range_generator<mtl::tag::major, compressed2D<t_int>>::type cType;
    typedef typename mtl::traits::range_generator<mtl::tag::nz,  cType>::type icType;
    typename mtl::traits::row<compressed2D<t_int>>::type row(K);
    typename mtl::traits::col<compressed2D<t_int>>::type col(K);
    typename mtl::traits::value<compressed2D<t_int>>::type value(K);
    for (cType cursor(mtl::begin<mtl::tag::major>(K)), cend(mtl::end<mtl::tag::major>(K)); cursor != cend; ++cursor )
        for (icType icursor(mtl::begin<mtl::tag::nz>(cursor)), icend(mtl::end<mtl::tag::nz>(cursor)); icursor != icend; ++icursor)
        {
            value(*icursor, c4n.size());
            c4n.emplace_back( (c4n[row(*icursor)] + c4n[col(*icursor)])/2 );
        }
    T = K + trans(K);
    swap(K,T);

    // Define new elements.
    t_n4e n4e_new; n4e_new.reserve( 4*n4e.size() );
    for ( t_int curElem = 0, nrElem = n4e.size(); curElem < nrElem; ++curElem )
    {
        n4e_new.push_back( std::array<t_int,3>{{n4e[curElem][0], K(n4e[curElem][0],n4e[curElem][1]), K(n4e[curElem][2],n4e[curElem][0])}} );
        n4e_new.push_back( std::array<t_int,3>{{K(n4e[curElem][0],n4e[curElem][1]), n4e[curElem][1], K(n4e[curElem][1],n4e[curElem][2])}} );
        n4e_new.push_back( std::array<t_int,3>{{K(n4e[curElem][0],n4e[curElem][2]), K(n4e[curElem][1],n4e[curElem][2]), n4e[curElem][2]}} );
        n4e_new.push_back( std::array<t_int,3>{{K(n4e[curElem][0],n4e[curElem][1]), K(n4e[curElem][1],n4e[curElem][2]), K(n4e[curElem][2],n4e[curElem][0])}} );
    }
    swap( n4e, n4e_new );

    // Define new Dirichlet boundary.
    t_n4ed dbnd_new; dbnd_new.reserve( 2*dbnd.size() );
    for ( t_int curEdge = 0, nrEdges = dbnd.size(); curEdge < nrEdges; ++curEdge )
    {
        dbnd_new.push_back( std::array<t_int,2>{{dbnd[curEdge][0], K(dbnd[curEdge][0],dbnd[curEdge][1])}} );
        dbnd_new.push_back( std::array<t_int,2>{{K(dbnd[curEdge][0],dbnd[curEdge][1]),dbnd[curEdge][1]}} );
    }
    swap( dbnd, dbnd_new );

    // Define new Neumann boundary.
    t_n4ed nbnd_new; nbnd_new.reserve( 2*nbnd.size() );
    for ( t_int curEdge = 0, nrEdges = nbnd.size(); curEdge < nrEdges; ++curEdge )
    {
        nbnd_new.push_back( std::array<t_int,2>{{nbnd[curEdge][0], K(nbnd[curEdge][0],nbnd[curEdge][1])}} );
        nbnd_new.push_back( std::array<t_int,2>{{K(nbnd[curEdge][0],nbnd[curEdge][1]),nbnd[curEdge][1]}} );
    }
    swap( nbnd, nbnd_new );
}

void comp_inv_det( t_matrix_3x3& P_inv, t_double& P_det, const t_matrix_3x3& P )
{
    P_det = P(0,0) * P(1,1) * P(2,2) +
            P(1,0) * P(2,1) * P(0,2) +
            P(0,1) * P(1,2) * P(2,0) -
            P(0,2) * P(1,1) * P(2,0) -
            P(0,1) * P(1,0) * P(2,2) -
            P(1,2) * P(2,1) * P(0,0);
    P_inv(0,0) = P(1,1)*P(2,2) - P(1,2)*P(2,1);
    P_inv(0,1) = P(0,2)*P(2,1) - P(0,1)*P(2,2);
    P_inv(0,2) = P(0,1)*P(1,2) - P(0,2)*P(1,1);
    P_inv(1,0) = P(1,2)*P(2,0) - P(1,0)*P(2,2);
    P_inv(1,1) = P(0,0)*P(2,2) - P(0,2)*P(2,0);
    P_inv(1,2) = P(0,2)*P(1,0) - P(0,0)*P(1,2);
    P_inv(2,0) = P(1,0)*P(2,1) - P(1,1)*P(2,0);
    P_inv(2,1) = P(0,1)*P(2,0) - P(0,0)*P(2,1);
    P_inv(2,2) = P(0,0)*P(1,1) - P(0,1)*P(1,0);
    P_inv /= P_det;
}

void fem15( const t_c4n& c4n, const t_n4e& n4e, const t_n4ed& dbnd, const t_n4ed&, t_u4n& u4n )
{
    const t_int nrNodes = c4n.size();
    compressed2D<t_double> A(nrNodes,nrNodes);
    dense_vector<t_double> b(nrNodes, 0.0), x(nrNodes);
    t_double P_det;
    t_matrix_3x3 P, P_inv, A_T;
    t_matrix_3x2 Q, T;
    set_to_zero(T); T(1,0) = 1.0; T(2,1) = 1.0;

    // Build the system of equations.
    {
        mat::inserter<compressed2D<t_double>,update_plus<t_double>> ins(A);
        for ( t_int curElem = 0, nrElem = n4e.size(); curElem < nrElem; ++curElem )
        {
            // Build local stiffness matrices A_T for piecewise linear shape functions and insert into the global stiffness matrix A.
            P(0,0) = 1.0; P(0,1) = 1.0; P(0,2) = 1.0;
            P(1,0) = c4n[n4e[curElem][0]](0); P(1,1) = c4n[n4e[curElem][1]](0); P(1,2) = c4n[n4e[curElem][2]](0);
            P(2,0) = c4n[n4e[curElem][0]](1); P(2,1) = c4n[n4e[curElem][1]](1); P(2,2) = c4n[n4e[curElem][2]](1);
            comp_inv_det(P_inv,P_det,P);
            Q = P_inv * T;
            A_T = (P_det / 2) * Q * trans(Q);
            std::vector<t_int> curDofIdx{ n4e[curElem][0], n4e[curElem][1], n4e[curElem][2] }; //  MTL4 does not support std::array yet.
            ins << element_matrix(A_T,curDofIdx);

            // Build global load vector for f = 1, g = 0, u_D = 0.
            for ( t_int m = 0; m < 3; ++m )
                b(n4e[curElem][m]) += P_det / 6;
        }
    }

    // Determine free nodes, i.e., exclude Dirichlet nodes.
    std::vector<bool> isFreeNode(nrNodes,true);
    for ( t_int i = 0, n = dbnd.size(); i < n; ++i )
        isFreeNode[dbnd[i][0]] = isFreeNode[dbnd[i][1]] = false;
    std::vector<t_int> freeNodes; freeNodes.reserve(nrNodes);
    for ( t_int i = 0; i < nrNodes; ++i )
        if ( isFreeNode[i] )
            freeNodes.push_back(i);
    mat::traits::permutation<>::type freeNodes_perm = mat::reorder( freeNodes, nrNodes );

    // Build reduced system of equations for free nodes.
    compressed2D<t_double> A_freeNodes( freeNodes_perm * A * trans(freeNodes_perm) );
    dense_vector<t_double> b_freeNodes( freeNodes_perm * b );
    dense_vector<t_double> x_freeNodes( num_cols(A_freeNodes), 0.0 );

    // Solve (reduced) problem equations Ax = b with left preconditioner PC
    pc::ic_0<compressed2D<t_double>> PC(A_freeNodes);
    basic_iteration<t_double> iter(b_freeNodes, 10000, 1.e-6);
    cg(A_freeNodes, x_freeNodes, b_freeNodes, PC, iter);
    x = trans(freeNodes_perm) * x_freeNodes;

    u4n.clear(); u4n.resize(nrNodes);
    std::copy( x.begin(), x.end(), u4n.begin() );
}

int main()
{
    // Square
    t_c4n c4n{ t_vector_2{0.0,0.0}, t_vector_2{1.0,0.0}, t_vector_2{1.0,1.0}, t_vector_2{0.0,1.0} };
    t_n4e n4e{ std::array<t_int,3>{{0,1,3}}, std::array<t_int,3>{{3,1,2}} };
    t_n4ed dbnd{ std::array<t_int,2>{{0,1}}, std::array<t_int,2>{{1,2}}, std::array<t_int,2>{{2,3}}, std::array<t_int,2>{{3,0}} };
    t_n4ed nbnd;

//    // Square 2
//    t_c4n c4n{ t_vector_2{-1.0,-1.0}, t_vector_2{0.0,-1.0}, t_vector_2{1.0,-1.0},
//               t_vector_2{-1.0,0.0}, t_vector_2{0.0,0.0}, t_vector_2{1.0,0.0},
//               t_vector_2{-1.0,1.0}, t_vector_2{0.0,1.0}, t_vector_2{1.0,1.0} };
//    t_n4e n4e{ std::array<t_int,3>{{0,1,4}}, std::array<t_int,3>{{0,4,3}},
//               std::array<t_int,3>{{1,2,4}}, std::array<t_int,3>{{4,2,5}},
//               std::array<t_int,3>{{4,5,8}}, std::array<t_int,3>{{4,8,7}},
//               std::array<t_int,3>{{3,4,6}}, std::array<t_int,3>{{6,4,7}} };
//    t_n4ed dbnd{ std::array<t_int,2>{{0,1}}, std::array<t_int,2>{{1,2}},
//                 std::array<t_int,2>{{2,5}}, std::array<t_int,2>{{5,8}},
//                 std::array<t_int,2>{{8,7}}, std::array<t_int,2>{{7,6}},
//                 std::array<t_int,2>{{6,3}}, std::array<t_int,2>{{3,0}} };
//    t_n4ed nbnd;

//    // L-shape
//    t_c4n c4n{ t_vector_2{-1.0,-1.0}, t_vector_2{0.0,-1.0},
//               t_vector_2{-1.0,0.0}, t_vector_2{0.0,0.0}, t_vector_2{1.0,0.0},
//               t_vector_2{-1.0,1.0}, t_vector_2{0.0,1.0}, t_vector_2{1.0,1.0} };
//    t_n4e n4e{ std::array<t_int,3>{{0,1,3}}, std::array<t_int,3>{{0,3,2}},
//               std::array<t_int,3>{{3,4,7}}, std::array<t_int,3>{{3,7,6}},
//               std::array<t_int,3>{{2,3,5}}, std::array<t_int,3>{{5,3,6}} };
//    t_n4ed dbnd{ std::array<t_int,2>{{0,1}}, std::array<t_int,2>{{1,3}},
//                 std::array<t_int,2>{{3,4}}, std::array<t_int,2>{{4,7}},
//                 std::array<t_int,2>{{7,6}}, std::array<t_int,2>{{6,5}},
//                 std::array<t_int,2>{{5,2}}, std::array<t_int,2>{{2,0}} };
//    t_n4ed nbnd;

    t_u4n u4n;
    for ( t_int curRef = 0, nrRef = 8; curRef < nrRef; ++curRef )
    {
        red_refine(c4n,n4e,dbnd,nbnd);
        std::cout << "nrNodes: " << c4n.size() << "\n";
        fem15(c4n,n4e,dbnd,nbnd,u4n);
        write_svg(c4n,n4e,u4n,"test" + std::to_string(curRef) + ".svg");
    }

    return 0;
}
