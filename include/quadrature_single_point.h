#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>
//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template <typename Ret, typename Integrand_Func>
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q,
                                    Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                    Integrand_Func integrand)
{
    // Select a point within the tetrahedron, we choose centroid
    Eigen::Vector3d centroid;
    for (int i = 0; i < 4; ++i)
    {
        centroid += q.segment<3>(3 * element(i));
    }
    centroid = centroid / 4;
    integrand(integrated, q, element, centroid);
    integrated = volume * integrated;


}
