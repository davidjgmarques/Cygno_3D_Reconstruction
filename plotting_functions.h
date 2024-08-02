#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void create_and_print_wf_graph_lines2 (std::string filename, std::vector<int> time, std::shared_ptr<std::vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2, std::vector<std::pair<int,double>> peaks);
void print_graph_lines2 (TGraph *graph, std::string title, std::string x_axis, std::string y_axis, double yMin, double yMax, TLine *l1, TLine *l2, TMarker *p1, TMarker *p2);

void addPoints_BAT_CAM( TCanvas *image, std::vector<std::pair<double,double>> points, std::string mode, std::string title);

void printTrackProfiles ( TH1D *h1, TH1D *h2, std::string title);

void addTracks( TCanvas *image, TH2F *histo, TH2F* track, int fminx, int fminy, std::string nometh2);

void build_3D_vector (double x0, double x1, double y0, double y1, double z0, double z1,
 double l_xy, double a_xy, double l_z, int d_z, double p_z, double a_z, double length,
 int ru, int pi, int tr);

#endif