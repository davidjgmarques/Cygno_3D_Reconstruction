
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>

#include "cl_trg_association.h"
#include "bat_functions.h"

void one_to_one_association(std::vector<std::tuple<double, size_t, size_t>> &dists, const std::vector<AlphaTrackPMT>& pmt_entries, const std::vector<AlphaTrackCAM>& cam_entries, bool verbose) {

    std::vector<std::pair<double,double>> sing_point_dists;
    double tot_avg_dist = 0;

    for (size_t i = 0; i < cam_entries.size(); ++i) {
        for (size_t j = 0; j < pmt_entries.size(); ++j) {

            tot_avg_dist = 0;
            sing_point_dists.clear();

            calculate_distance(pmt_entries[j].track_pmt, cam_entries[i].track_cam, sing_point_dists, verbose);

            for ( const auto point : sing_point_dists ) {

                // Sum of *all* distances is X and Y
                // The points are alreadu ordered by X
                tot_avg_dist += std::sqrt( point.first * point.first + point.second * point.second);
            } 

            tot_avg_dist /= sing_point_dists.size();

            dists.emplace_back(tot_avg_dist, i, j);
        }
    }

    std::sort(dists.begin(), dists.end());
}