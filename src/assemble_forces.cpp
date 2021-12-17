#include <assemble_forces.h>
#include <iostream>
#include <dV_linear_tetrahedron_dq.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  f - the vector 3n vector of forces acting on each node of the mass-spring system

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D, const std::string&  energy_type) {

    f.resize(3  * q.rows());
    f.setZero();

    // Iterate over all tetrahedrons in the connectivity matrix
    for (int i = 0; i < T.rows(); ++i)
    {
        Eigen::Vector12d dV_dq;
        Eigen::RowVector4i element = T.row(i);
        double volume = v0(i);
        dV_linear_tetrahedron_dq(dV_dq, q, V, element, volume, C, D, energy_type);

        // Compute force between every spring pair within the tetrahedron
        for (int j = 0; j < 4; j++) {
            f(3 * T(i, j))   -= dV_dq(3 * j);
            f(3 * T(i, j)+1) -= dV_dq(3 * j+1);
            f(3 * T(i, j)+2) -= dV_dq(3 * j+2);
        }
    }

    };