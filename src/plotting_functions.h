#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void create_and_print_wf_graph_lines3 (std::string filename, std::vector<int> time, std::shared_ptr<std::vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2,  std::vector<std::pair<int,double>> peaks_crown,  std::vector<std::pair<int,double>> peaks_energy_dep);
void print_graph_lines3 (TGraph *graph, std::string title, std::string x_axis, std::string y_axis, double yMin, double yMax, TLine *l1, TLine *l2, std::vector<TMarker*> markers_cr, std::vector<TMarker*> markers_ed);

void Points_BAT_CAM(std::vector<std::pair<double,double>> batPoints, std::vector<std::pair<double,double>> camPoints, std::string title);

void printTrackProfilesAndFit ( TH1D *h1, TH1D *h2, std::string title, TFitResultPtr &fitResult, bool save);

void doubleGaussianTest( TH1D *h1, TH1D *h2, std::string title, TFitResultPtr &fitResult, bool save);

void centralGaussianTest( TH1D *h1, TH1D *h2, std::string title, TFitResultPtr &fitResult, bool save);

double getHistogramRMS(TH1D* histo);

double getTrackProfileWidth(TH1D* histo, double percentage, int consecutiveBins,bool verb);

void addTracks( TCanvas *image, TH2F *histo, TH2F* track, int fminx, int fminy, std::string nometh2);

void build_3D_vector (double x0, double x1, double y0, double y1, double z0, double z1,
 double l_xy, double a_xy, double l_z, int d_z, double p_z, double a_z, double length,
 int ru, int pi, int tr, bool cloud, double sigma);

#endif