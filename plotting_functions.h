#ifndef FUNCTIONS_H
#define FUNCTIONS_H

using namespace std;

void create_and_print_wf_graph_lines2 (string filename, vector<int> time, shared_ptr<vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2, vector<pair<int,double>> peaks);
void print_graph_lines2 (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax, TLine *l1, TLine *l2, TMarker *p1, TMarker *p2);


#endif