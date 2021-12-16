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
//    psi = C * (std::pow(J, -2. / 3.) * (F.transpose() * F).trace() - 3) + D * std::pow(J - 1, 2);

    // Stable Neo-hookean
    double mu = C * 2;
    double lambda = D * 2;
    double alpha = 1 + mu / lambda  - mu / (4 * lambda);
    psi = mu * 0.5 * ((F.transpose() * F).trace() - 3) + lambda * 0.5 * std :: pow(J - alpha, 2) - mu * 0.5 * log((F.transpose() * F).trace() + 1);



}