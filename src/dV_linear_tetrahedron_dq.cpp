#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

//Input:
//  q - generalized coordinates for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  v0 - the undeformed tetrahedron volume
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  dV - the 12x1 gradient of the potential energy for a single tetrahedron

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

//   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q,
//   Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
//
//    };
//    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);

    auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

        // Calculate F
        Eigen::Matrix3d F;
        Eigen::Matrix34d t;
        for (int i = 0; i < 4; ++i)
        {
            t.col(i) = q.segment<3>(3 * element(i));
        }
        Eigen::Matrix43d dphi;
        dphi_linear_tetrahedron_dX(dphi, V, element, X);
        F = t * dphi;

        // Compute dPsi/dF
        Eigen::Vector9d dpsi_by_dF;
        dpsi_neo_hookean_dF(dpsi_by_dF, F, C, D);

        // Build B of F = B * x
        Eigen::Matrix<double, 9, 12> B;
        B.setZero();
        // Iteratively set 1x3 column vectors in B using row vectors in dphi
        for (int i = 0; i < 4; ++i)
        {
            B.block(0, 0 + 3 * i, 3, 1) = dphi.row(i).transpose();
            B.block(3, 1 + 3 * i, 3, 1) = dphi.row(i).transpose();
            B.block(6, 2 + 3 * i, 3, 1) = dphi.row(i).transpose();
        }
        dV = B.transpose() * dpsi_by_dF;
    };

    // This function returns vol * lambda(x), where x is a point within the tetrahedron (we choose centroid)
    // Here, lambda(x) = E.T * B.T * dpsi_dF
    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);


}