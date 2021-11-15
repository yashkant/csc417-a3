#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>
typedef Eigen::Triplet<double> Trip;

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  V_skin - lx3 matrix of vertices of the display mesh
//Output:
//  N - the lxn sparse skinning matrix

bool same_side(Eigen::Vector3d &v1, Eigen::Vector3d &v2, Eigen::Vector3d &v3,
        Eigen::Vector3d &v4, Eigen::Vector3d &p){
            Eigen::Vector3d normal = (v2 - v1).cross(v3-v1);
            double dot_v4 = normal.dot(v4-v1);
            double dot_p = normal.dot(p - v1);
            bool same_side = dot_p * dot_v4 > 0;
            return same_side;

        }

bool in_tet(Eigen::Vector3d &v1, Eigen::Vector3d &v2, Eigen::Vector3d &v3,
        Eigen::Vector3d &v4, Eigen::Vector3d &p){
            bool ss1 = same_side(v1, v2, v3, v4, p);
            bool ss2 = same_side(v2, v3, v4, v1, p);
            bool ss3 = same_side(v3, v4, v1, v2, p);
            bool ss4 = same_side(v4, v1, v2, v3, p);
            bool in_tet = ss1 && ss2 && ss3 && ss4;
            return in_tet;
        }


void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V,
                            Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::MatrixXd> V_skin)
{
    std::vector<Trip> triplets;

    // For every skin vertex
    for (int i = 0; i < V_skin.rows(); ++i){
        Eigen::Vector3d sv = V_skin.row(i);
        // For every tetrahedron
        for (int j = 0; j < T.rows(); ++j){
            Eigen::Vector3d v1 = V.row(T(j, 0));
            Eigen::Vector3d v2 = V.row(T(j, 1));
            Eigen::Vector3d v3 = V.row(T(j, 2));
            Eigen::Vector3d v4 = V.row(T(j, 3));
            bool in_tetr = in_tet(v1, v2, v3, v4, sv);
            if (in_tetr) {
                Eigen::Vector4d phi;
                phi_linear_tetrahedron(phi, V, T.row(j), V_skin.row(i).transpose());
                for (int k = 0; k < 4; ++k)
                    triplets.push_back(Trip(i, T(j, k), phi(k)));
            };
        };
    }
    std::cout<<"Build skinning matrix!"<<"\n";
    N.resize(V_skin.rows(), V.rows());
    N.setFromTriplets(triplets.begin(), triplets.end());


}