#include <dV_spring_particle_particle_dq.h>

//Input:
//  q0 - the generalized coordinates of the first node of the spring
//  q1 - the generalized coordinates of the second node of the spring
//  l0 - the undeformed length of the spring
//  stiffness - the stiffness constant for this spring
//Output:
//  f - the 6x1 per spring generalized force vector

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,
    Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {

    // Calculate deformation
    double deformation;
    deformation = (q1 - q0).norm() - l0;

    // Force magnitude
    // Force acts along the line joining q1 and q0
    double f_mag = stiffness * deformation;

    // Force on particle q0, will be vector going towards q1 from q0
    Eigen::Vector3d f_vec;
    f_vec.resize(q0.size());
    f_vec = -1 * (q1 - q0) / (q1-q0).norm();
    f.segment<3>(0) = f_vec * f_mag;

    // Force on particle q1, will be vector going towards q0 from q1
    f_vec = -1 * (q0 - q1) / (q0-q1).norm();
    f.segment<3>(3) = f_vec * f_mag;

    
}