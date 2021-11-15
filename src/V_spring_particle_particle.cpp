#include <V_spring_particle_particle.h>

//Input:
//  q0 - the generalized coordinates of the first node of the spring
//  q1 - the generalized coordinates of the second node of the spring
//  l0 - the undeformed length of the spring
//  stiffness - the stiffness constant for this spring
//Output:
//  V - potential energy for the spring


//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    // Calculate deformation
    double deformation;
    deformation = (q1 - q0).norm() - l0;
    // Calculate PE
    // V(q) = 0.5 * k * (sqrt(q^T x B^T x B x q) - l0)^2
    V = 0.5 * stiffness * std :: pow(deformation,2);

}