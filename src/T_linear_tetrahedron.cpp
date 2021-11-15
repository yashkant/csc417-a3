#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

//Input:
// qdot - generalied velocity of FEM system
// element - vertex indices of the element
// density - material density
// volume - volume of tetrahedron
//Output:
//  T - kinetic energy of tetrahedron

void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot,
 Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {

    // Kinetic energy of a tetrahedron is given by 0.5 * v * M * v.transpose()
    // Where M is the per tetrahedron mass-matrix and v is a 12-D vector
    Eigen::Matrix1212d M_per_tetrahedron;
    mass_matrix_linear_tetrahedron(M_per_tetrahedron, qdot, element, density, volume);

    Eigen::VectorXd qdot_tetrahedron;
    qdot_tetrahedron.resize(12);

    // Populate qdot_tetrahedron by selection
    for (int i=0; i < 4; i++){
        qdot_tetrahedron.segment(i, 3) = qdot.segment(element(i) * 3, 3);
    }

    // Build KE
    T = 0.5 * qdot_tetrahedron.transpose() * M_per_tetrahedron * qdot_tetrahedron;
}
