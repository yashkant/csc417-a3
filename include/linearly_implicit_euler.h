#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template <typename FORCE, typename STIFFNESS>
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt,
                                    const Eigen::SparseMatrixd &mass, FORCE &force, STIFFNESS &stiffness,
                                    Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness)
{

    // We are trying to solve (M-dt*dt*K) qdot_t+1 = M * qdot_t + dt * force(t)
    // This gives us an equation of form Ax = v, where x = qdot_t+1, and v = RHS, and A = (M-dt*dt*K)
    // Instead of taking inverse of A, we use LDLT
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;

    // Calculating RHS -- v = M * qdot_t + dt * force(t)
    force(tmp_force, q, qdot);
    Eigen::VectorXd v = mass * qdot + dt * tmp_force;

    // Calculating LHS -- A = (M-dt*dt*K) and x = qdot_t+1
    stiffness(tmp_stiffness, q, qdot);
    Eigen::SparseMatrixd A = mass - dt * dt * tmp_stiffness;

    // Solve by ldlt
    ldlt.compute(A);
    qdot = ldlt.solve(v);

    // Compute q from qdot
    q = q + dt * qdot;


}