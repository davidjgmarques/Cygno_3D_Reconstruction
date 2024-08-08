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

    double IP_X_cm;
    double IP_Y_cm;

    bool its_alpha;

    std::vector<std::pair<double,double>> track_cam;

    double energy;
    double nhits;
    double width;
    double xmean;
    double ymean;
    double rms;
    double tgausssigma;

};

struct AlphaTrackPMT {

    int run;
    int pic;
    int trg;

    int dir;        // -1 = towards GEM ; 1 = towards cathode; 0 = ambiguous
    double prob;    
    double trv_Z;
    int quad;

    bool its_alpha;

    std::vector<std::pair<double,double>> track_pmt;

    double energy;

    int num_peaks;
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

#endif // CL_TRG_ASSOCIATION_H