#include <fixed_point_constraints.h>
#include <algorithm>
typedef Eigen::Triplet<double> Trip;

//Input:
//  q_size - total number of scalar generalized coordinates (3n for n vertices)
//  indices - m indices (row ids in V) for fixed vertices
//Output:
//  P - 3(n-m)x3n sparse matrix which projects out fixed vertices

void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices)
{

    // Resize to (3(n-m), 3n)
    P.resize(q_size - 3 * indices.size(), q_size);

    // Weed out indices of fixed points and fill in the rest with 1s
    unsigned int num_inds;
    num_inds = q_size / 3;

    int row_count;
    row_count = 0;
    std::vector<Trip> trips;
    trips.reserve(3);
    for (int i = 0; i < num_inds; i += 1){
        // if this index is not fixed add corresponding entries to P
        if (std::find(indices.begin(), indices.end(), i) == indices.end()) {

            // Add selection values
            trips.push_back(Trip(row_count * 3, i * 3, 1));
            trips.push_back(Trip(row_count * 3 + 1, i * 3 + 1, 1));
            trips.push_back(Trip(row_count * 3 + 2, i * 3 + 2, 1));

            // Assign values
            P.setFromTriplets(trips.begin(), trips.end());

            // Move to next row
            row_count = row_count + 1;
        }
    }
}