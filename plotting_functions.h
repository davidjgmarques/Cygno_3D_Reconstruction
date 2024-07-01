#ifndef FUNCTIONS_H
#define FUNCTIONS_H

using namespace std;

void create_and_print_wf_graph_lines3 (string filename, vector<int> time, shared_ptr<vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2, double p1t, double p1a, double p2t, double p2a);
void print_graph_lines3 (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax, TLine *l1, TLine *l2, TMarker *p1, TMarker *p2);


#endif