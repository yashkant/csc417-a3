#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template<typename Ret, typename Integrand_Func>

// Usage in V_linear_tetrahedron:
//      quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);

inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q, 
                                               Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                               Integrand_Func integrand) {
    // Select a point within the tetrahedron
    Eigen::Vector3d centroid;
    for (int i = 0; i < 4; ++i)
    {
        centroid += q.segment(3 * element(i), 3);
    }
    centroid = centroid / 4;
    integrand(integrated, q, element, centroid);
    integrated = volume * integrated;


}

