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
                     double C, double D, const std::string& energy_type) {

    double J = F.determinant();
    double mu = C * 2;
    double lambda = D * 2;
    double alpha = 1 + mu / lambda  - mu / (4 * lambda);

    if (energy_type == "smith_14"){
        // Stable Neo-hookean (Smith and Kim, 2018, Equation 14, https://dl.acm.org/doi/10.1145/3180491)
        psi = mu * 0.5 * ((F.transpose() * F).trace() - 3) + lambda * 0.5 * std :: pow(J - alpha, 2) - mu * 0.5 * log((F.transpose() * F).trace() + 1);
    }
    else if (energy_type == "smith_13"){
        // Stable Neo-hookean (Smith and Kim, 2018, Equation 13, https://dl.acm.org/doi/10.1145/3180491)
        // Without meta-stability
        psi = mu * 0.5 * ((F.transpose() * F).trace() - 3) + lambda * 0.5 * std :: pow(J - 1, 2) - mu * (J-1);
    }
    else if (energy_type == "bower"){
        // Neo-hookean strain energy (Bower, 2009)
        psi = mu * 0.5 * (std::pow(J, -2. / 3.) * (F.transpose() * F).trace() - 3) + lambda * 0.5 * std::pow(J - 1, 2);
    }
    else if (energy_type == "wang"){
        // Neo-hookean strain energy (Wang and Yang, 2016)
        psi = mu * 0.5 * (std::pow(J, -2. / 3.) * (F.transpose() * F).trace() - 3) + lambda * 0.5 * std::pow(J - 1, 1);

    }
    else if (energy_type == "ogden"){
        // Neo-hookean strain energy (Ogden, 1984)
        psi = mu * 0.5 * ((F.transpose() * F).trace() - 3) - mu * log(J) + 0.5 * lambda * std::pow(log(J), 2);
    }
    else{
        exit(0);
    }

}