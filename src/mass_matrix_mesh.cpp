#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

typedef Eigen::Triplet<double> Trip;

//Input:
//  qdot - generalized velocity for the FEM system
//  T - the mx4 vertex indices for tet mesh
//  density - density of material
//  v0 - the undeformed tetrahedra volumes
//Output:
//  M - Sparse mass matrix for the whole mesh.

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot,
Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

    // Resize the Mass Matrix M to 3N x 3N
    M.resize(qdot.rows(), qdot.rows());

    // We assemble the Mass Matrix M from individual tetrahedron mass matrices
    // Vector to store row-col to value map for the new sparse Mass Matrix M
    std::vector<Trip> tetrahedron_triplets;
    tetrahedron_triplets.reserve(T.rows() * 12 * 12);

    // Iterate over every single tetrahedron
    for (int tet = 0; tet < T.rows(); ++tet)
    {
        // Build per-tetrahedron mass matrix
        Eigen::Matrix1212d Mtet;
        Eigen::RowVectorXi element = T.row(tet);
        double volume = v0(tet);
        mass_matrix_linear_tetrahedron(Mtet, qdot, element, density, volume);

        // Iterate over every single pair of vertices of the tetrahedron
        for (int v1 = 0; v1 < 4; ++v1)
        {
            for (int v2 = 0; v2 < 4; ++v2)
            {
                // Iterate over every single element of corresponding mass matrix
                // between v1 and v2
                for (int row = 0; row < 3; ++row)
                {
                    for (int col = 0; col < 3; ++col)
                    {
                        tetrahedron_triplets.push_back(Trip(3 * element(v1) + row, 3 * element(v2) + col, Mtet(3 * v1 + row, 3 * v2 + col)));
                    }
                }
            }
        }
    }
    M.setFromTriplets(tetrahedron_triplets.begin(), tetrahedron_triplets.end());
}