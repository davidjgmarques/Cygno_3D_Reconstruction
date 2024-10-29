#ifndef CL_TRG_ASSOCIATION_H
#define CL_TRG_ASSOCIATION_H

#include <vector>
#include <tuple>

struct AlphaTrackCAM {

    int run;
    int pic;
    int cluster;

    double angle_XY; 
    double trv_XY;
    int quad;

    double begin_X_cm;
    double begin_Y_cm;

    double end_X_cm;
    double end_Y_cm;

    bool its_alpha;

    std::vector<std::pair<double,double>> track_cam;

    double energy;
    double nhits;
    double width;
    double xmean;
    double ymean;
    double rms;
    double tgausssigma;

    double fitSig;
    double calc_abs_Z;

};

struct AlphaTrackPMT {

    int run;
    int pic;
    int trg;

    int dir;        // -1 = towards cathode; 1 = towards GEM; 0 = ambiguous
    double prob;    
    double trv_Z;
    int quad;

    bool its_alpha;

    std::vector<std::pair<double,double>> track_pmt;

    double energy;

    double num_peaks;
};

// This structure is later useful to group the alpha by {run,event}
struct RunPicKey {
    int run;
    int pic;

    bool operator<(const RunPicKey& other) const {
        return std::tie(run, pic) < std::tie(other.run, other.pic);
    }
};

void one_to_one_association(std::vector<std::tuple<double, size_t, size_t>> &dists, const std::vector<AlphaTrackPMT>& pmt_entries, const std::vector<AlphaTrackCAM>& cam_entries, bool verbose);

bool found_clusters_in_evt(const std::vector<AlphaTrackCAM>& CAM_alphas, int pmt_run, int pmt_event);

void deleteNonAlphaDirectories(const char* filename, bool deleteAll);

void angle_3D_reverse(double angle_cam, std::vector<std::pair<double, double>> &points_cam);

double estimate_absolute_Z(double sigma);

int generate_random_direction();



#endif // CL_TRG_ASSOCIATION_H