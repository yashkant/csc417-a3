 
 #include <mass_matrix_linear_tetrahedron.h>

//Input:
//  qdot - generalized velocity for the FEM system
//  element - the 1x4 vertex indices for this tetrahedron
//  density - density of material
//  volume - the undeformed tetrahedron volume
//Output:
//  M - dense 12x12 per-tetrahedron mass matrix

 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot,
 Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {

    // Assemble a 12x12 Mass Matrix, by solving Integration (Phi_r(X) * Phi_s(X) * d_omega * I)
    // Turns out that M is a symmetric matrix, with same elements on diagonal and off-diagonal
    // respectively. Diagonal element are -- 6 * rho * vol * (1/60) and Off-diagonal are --
    // 6 * rho * vol * (1/120)

    Eigen::MatrixXd I;
    I.resize(3,3);
    I.setIdentity();

    auto mass = density * volume;
    auto diag = (1.0 / 10.0) * mass;
    auto off_diag = (1.0 / 20.0) * mass;

    for (int row=0; row<4; ++row){
        for (int col=0; col<4; ++col){
            // Diagonal
            if (row == col){
                M.block(row * 3, col * 3, 3, 3) = diag * I;
            }
            // Off-diagonal
            else {
                M.block(row * 3, col * 3, 3, 3) = off_diag * I;
            }
        }
    }
 }