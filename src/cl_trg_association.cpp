
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <TFile.h>
#include <TKey.h>

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

bool found_clusters_in_evt(const std::vector<AlphaTrackCAM>& CAM_alphas, int pmt_run, int pmt_event) {

    bool found = false;
    for (const auto &cam_alpha : CAM_alphas) {

        if (cam_alpha.run == pmt_run && cam_alpha.pic == pmt_event) {
            found = true;
            break;
        }
    }

    return found;    
}

void angle_3D_reverse(double angle_cam, std::vector<std::pair<double, double>> &points_cam)
{
    // Check the angle and reverse the vector if necessary
    // Only relevant for plotting purposes, since the angle is saved directly in the structure
    // Irrelevant for BAT-CAM distance matching, because there  I actually order the points by X.
    // Added later due to the change in the way the track end point is calculated.

    if (angle_cam > 90.0 || angle_cam < -90.0)
    {
        std::reverse(points_cam.begin(), points_cam.end());
    }
}

void deleteNonAlphaDirectories(const char* filename, bool deleteAll) {

    std::cout << "Deleting non-alpha events from root file..." << std::endl;

    // Open the ROOT file
    TFile* file_root = TFile::Open(filename, "UPDATE");
    if (!file_root || file_root->IsZombie()) {
        std::cout << "Error opening file: " << filename << std::endl;
        return;
    }

    // Get the list of directories
    TIter next(file_root->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        // Check if the key is a directory
        if (strcmp(key->GetClassName(), "TDirectoryFile") != 0) continue;

        // Get the directory
        TDirectory* dir = (TDirectory*)file_root->Get(key->GetName());
        if (!dir) continue;

        // If deleteAll is true, delete the directory
        if (deleteAll) {
            file_root->cd();
            file_root->rmdir(dir->GetName());
            continue;
        } else {

            // Check if the directory contains "3D_vector"
            bool contains3DVector = false;
            TIter nextObj(dir->GetListOfKeys());
            TKey* objKey;
            while ((objKey = (TKey*)nextObj())) {
                std::string objName = objKey->GetName();
                if (objName.find("3D_vector") != std::string::npos) {
                    contains3DVector = true;
                    break;
                }
            }

            // If the directory doesn't contain "3D_vector", delete it
            if (!contains3DVector) {
                file_root->cd();
                file_root->rmdir(dir->GetName());
            }
        }
    }
}