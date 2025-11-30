/* -----------------------------------------------------------------------------
 - cl_trg_association.cpp
 - - Association utilities for matching camera (CAM) clusters with PMT (BAT)
 - - reconstructions. Functions here compute distances, perform one-to-one
 - - association ranking, provide simple geometric helpers and file-cleaning
 - - utilities used by the main analysis pipeline. No behavior changes.
 - ----------------------------------------------------------------------------- */

#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <TFile.h>
#include <TKey.h>
#include <random>

#include "cl_trg_association.h"
#include "bat_functions.h"

/* -----------------------------------------------------------------------------
 - one_to_one_association
 - - Compute average pointwise Euclidean distance between each CAM cluster
 - - and each PMT/BAT track and append tuples (avg_distance, cam_index, pmt_index)
 - - to the provided 'dists' vector. The output vector is sorted ascending by
 - - distance so the smallest-distance associations are first.
 - - Parameters:
 - -   dists: output vector of tuples (distance, cam_idx, pmt_idx)
 - -   pmt_entries: list of PMT/BAT reconstructed alpha tracks
 - -   cam_entries: list of camera reconstructed alpha tracks
 - -   verbose: enable per-point debug printing passed to calculate_distance
 - ----------------------------------------------------------------------------- */
void one_to_one_association(std::vector<std::tuple<double, size_t, size_t>> &dists, const std::vector<AlphaTrackPMT>& pmt_entries, const std::vector<AlphaTrackCAM>& cam_entries, bool verbose) {

    std::vector<std::pair<double,double>> dists_CAM_BAT;
    double tot_avg_dist = 0;

    for (size_t i = 0; i < cam_entries.size(); ++i) {
        for (size_t j = 0; j < pmt_entries.size(); ++j) {

            tot_avg_dist = 0;
            dists_CAM_BAT.clear();

            calculate_distance(pmt_entries[j].track_pmt, cam_entries[i].track_cam, dists_CAM_BAT, verbose);

            for ( const auto dist : dists_CAM_BAT ) {

                // Sum of *all* distances is X and Y
                // The dists had already been ordered by X
                // dist.first = (x2 - x1); dist.second = (y2 - y1)
                
                tot_avg_dist += std::sqrt( dist.first * dist.first + dist.second * dist.second);
            } 

            tot_avg_dist /= dists_CAM_BAT.size();

            dists.emplace_back(tot_avg_dist, i, j);
        }
    }

    std::sort(dists.begin(), dists.end());
}

/* -----------------------------------------------------------------------------
 - found_clusters_in_evt
 - - Check if any CAM alpha entries correspond to the given run/picture
 - - (used to decide whether to run PMT analysis for that event).
 - - Returns true if at least one CAM alpha has matching run and picture.
 - ----------------------------------------------------------------------------- */
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

/* -----------------------------------------------------------------------------
 - angle_3D_reverse
 - - Reverse the order of track points when the camera angle indicates the
 - - vector should be flipped for plotting. This does not change analysis
 - - associations (which rely on X-ordered points) and is purely cosmetic.
 - ----------------------------------------------------------------------------- */
void angle_3D_reverse(double angle_cam, std::vector<std::pair<double, double>> &points_cam)
{

    if (angle_cam > 90.0 || angle_cam < -90.0)
    {
        std::reverse(points_cam.begin(), points_cam.end());
    }
}

/* -----------------------------------------------------------------------------
 - estimate_absolute_Z
 - - Estimate absolute Z position from the measured transverse sigma using
 - - a diffusion model: Z = (sigma^2 - sigma0^2) / sigma_transverse^2.
 - - Returns -1 on invalid sigma or when the estimate is unphysical.
 - ----------------------------------------------------------------------------- */
double estimate_absolute_Z(double sigma) {

    // Z can be estimated using the track diffusion, from:
    // Z = (sigma_measured^2 - sigma_0^2) / sigma_transverse^2

    // From Lemon paper -- 500 V/cm
    // double sigma_zero       = 292.0E-4;             //um->cm. cm/sqrt(cm)  
    // double sigma_transverse = 130.0E-4;             //um->cm. cm/sqrt(cm)

    // From simulation -- Most up-to-date values
    // s_0 = 350 um/sqrt(cm)
    // s_t = 115 um/sqrt(cm) 

    // From MANGO
    // double sigma_zero = 390.0E-4;             //um->cm. cm/sqrt(cm) //from MANGO
    
    // From Rita Roque's PhD thesis
    // double sigma_zero = 500.0E-4;             //um->cm. cm/sqrt(cm)  
    // double sigma_transverse = 110.0E-4;       //um->cm. cm/sqrt(cm)

    // From a posteriori correction
    // double sigma_zero = 900.0E-4;             //um->cm. cm/sqrt(cm)  
    // double sigma_transverse = 115.0E-4;       //um->cm. cm/sqrt(cm)

    // Final from MANGO
    double sigma_zero = 780.0E-4;             //um->cm. cm/sqrt(cm)  
    double sigma_transverse = 130.0E-4;       //um->cm. cm/sqrt(cm)

    double Z = 0;

    if (sigma == 0) Z = -1;
    else            Z = (pow(sigma,2) - pow(sigma_zero,2)) / pow(sigma_transverse,2); 

    if (Z > 1E3) Z = -1;        //many fits don't work

    return Z;    
}

/* -----------------------------------------------------------------------------
 - deleteNonAlphaDirectories
 - - Open the output ROOT file and remove directories that do not contain
 - - a 3D_vector object (used to prune non-alpha output folders). If
 - - deleteAll is true, all directories are removed.
 - - Note: operates in-place on the provided ROOT file name.
 - ----------------------------------------------------------------------------- */
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

/* -----------------------------------------------------------------------------
 - generate_random_direction
 - - Return either -1 or +1 with equal probability. Used when head-tail
 - - determination is ambiguous and a random choice is preferable to losing
 - - statistics for direction-insensitive quantities.
 - ----------------------------------------------------------------------------- */
int generate_random_direction() {

    std::random_device rd;  // Seed
    std::mt19937 gen(rd()); // Random number generator
    std::uniform_int_distribution<> dis(0, 1); // Distribution that generates 0 or 1
    int random_dir = dis(gen) * 2 - 1; // Convert 0 or 1 to -1 or 1

    return random_dir;
}