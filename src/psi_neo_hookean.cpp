#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>


//Input:
//  F - the dense 3x3 deformation gradient
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  psi - the neohookean energy

void psi_neo_hookean(double &psi,
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D) {

    // Neo-hookean strain energy density as mentioned in the readme
    double J = F.determinant();
    psi = C * (std::pow(J, -2. / 3.) * (F.transpose() * F).trace() - 3) + D * std::pow(J - 1, 2);


}