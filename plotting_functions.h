#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void create_and_print_wf_graph_lines2 (std::string filename, std::vector<int> time, std::shared_ptr<std::vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2, std::vector<std::pair<int,double>> peaks);
void print_graph_lines2 (TGraph *graph, std::string title, std::string x_axis, std::string y_axis, double yMin, double yMax, TLine *l1, TLine *l2, TMarker *p1, TMarker *p2);
void print_BAT_CAM_match( TH2F *h_track, std::vector<std::pair<double,double>> points, std::string mode, std::string title);
void printTrackProfiles ( TH1D *h1, TH1D *h2, std::string title);

#endif