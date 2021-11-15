#include <assemble_stiffness.h>
typedef Eigen::Triplet<double> Trip;

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  K - the sparse, global stiffness matrix

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::VectorXd> qdot,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
        double C, double D) {

    K.resize(q.rows(), q.rows());
    std::vector<Trip> trips;
    trips.reserve(T.rows() * 4 * 4 * 9);

    // Iterate over all tetrahedrons in the connectivity matrix
    for (int i = 0; i < T.rows(); ++i)
    {
        Eigen::Matrix1212d d2V_dq2;

        // Calculate per tetrahedron stiffness matrix
        d2V_linear_tetrahedron_dq2(d2V_dq2, q, V, T.row(i), v0(i), C, D);

        // Fit negative hessian values for every pair of tetrahedron vertices in K
        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < 4; ++k)
            {
                for (int row = 0; row < 3; ++row)
                {
                    for (int col = 0; col < 3; ++col)
                    {
                        trips.push_back(Trip(3 * T(i, j) + row, 3 * T(i, k) + col, -d2V_dq2(3 * j + row, 3 * k + col)));
                    }
                }
            }
        }
    }
    K.setFromTriplets(trips.begin(), trips.end());
    };
