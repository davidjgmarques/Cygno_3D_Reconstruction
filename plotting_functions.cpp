#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TLine.h"
#include "TAxis.h"

#include "plotting_functions.h"
using namespace std;


void create_and_print_wf_graph_lines2 (string filename, vector<int> time, shared_ptr<vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2, vector<pair<int,double>> peaks) {

    TGraph *gWaveform = new TGraph();
    string newname = filename + "";

    TLine * line1 = new TLine(start1,level1,end1,level1);
    TLine * line2 = new TLine(start2,level2,end2,level2);

    TMarker *marker1 = new TMarker(peaks[0].first, peaks[0].second, 43); // Point (6, 36) with style 20
    TMarker *marker2 = new TMarker(peaks[1].first, peaks[1].second, 43); // Point (6, 36) with style 20

    for (int k = 0; k < time.size(); k++){

        gWaveform -> SetPoint ( k, time[k], (*ampl)[k]);
    }
    print_graph_lines2(gWaveform, newname, "Sample [#]", "ADC counts [#]", 0, 4000, line1, line2, marker1, marker2);
    
    // I need to save the canvas to also print the lines
    // gWaveform->SetName(newname.c_str());
    // gWaveform->Write(newname.c_str(),TObject::kWriteDelete);
}

void print_graph_lines2 (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax, TLine *l1, TLine *l2, TMarker *p1, TMarker *p2){

    TCanvas *c = new TCanvas("","", 800, 600);
    c->cd();
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(x_axis.c_str());
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetXaxis()->SetTitleOffset(1);
    graph->GetYaxis()->SetTitle(y_axis.c_str());
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleOffset(1);
    graph->SetLineColor(kAzure-5);
    graph->SetLineWidth(3);
    graph->SetMarkerColor(kAzure-5);
    // graph->GetYaxis()->SetRangeUser(yMin,yMax);
    graph->Draw("apl");

    l1->SetLineColor(kRed-7);
    l1->SetLineWidth(3);
    l1->SetLineStyle(9);
    l1->Draw("same");

    l2->SetLineColor(kGray+1);
    l2->SetLineWidth(3);
    l2->SetLineStyle(9);
    l2->Draw("same");

    p1->SetMarkerColor(kOrange); 
    p2->SetMarkerColor(kOrange); 
    p1->SetMarkerStyle(23);
    p2->SetMarkerStyle(23);
    p1->SetMarkerSize(1.5);
    p2->SetMarkerSize(1.5);
    p1->Draw("same");
    p2->Draw("same");

    c->SetName(title.c_str());
    c->Write(title.c_str(),TObject::kWriteDelete);
}