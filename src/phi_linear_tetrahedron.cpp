#include <phi_linear_tetrahedron.h>

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the underformed space at which to compute the energy density
//Output:
//  phi - the 4x1 values the basis functions

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V,
Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {

    // Extract coordinates of tetrahedron corners
    auto X0 = V.row(element(0)).transpose();
    auto X1 = V.row(element(1)).transpose();
    auto X2 = V.row(element(2)).transpose();
    auto X3 = V.row(element(3)).transpose();

    // Build Matrix T
    Eigen::Matrix3d T;
    T.col(0) = X1 - X0;
    T.col(1) = X2 - X0;
    T.col(2) = X3 - X0;

    // Phi(1,2,3) = T^inv * (X-X0)
    phi.segment(1, 3) = T.inverse() * (x-X0);
    phi(0) = 1 - phi(1) - phi(2) - phi(3);

}