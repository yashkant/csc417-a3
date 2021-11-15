#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

//Input:
// q - generalized coordinates of FEM system
// V - vertex matrix for the mesh
// element - vertex indices of the element
// volume - volume of tetrahedron
// C,D - material parameters
//Output:
//  energy - potential energy of tetrahedron


void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

                                                                                                             auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd> q,
     Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        Eigen::Matrix3d F;
        Eigen::Matrix34d t;
        for (int i = 0; i < 4; ++i)
        {
            t.col(i) = q.segment<3>(3 * element(i));
        }
        Eigen::Matrix43d dphi;
        dphi_linear_tetrahedron_dX(dphi, V, element, X);
        F = t * dphi;

        psi_neo_hookean(e, F, C, D);
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);

}