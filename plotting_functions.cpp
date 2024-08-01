#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TLine.h"
#include "TAxis.h"
#include "TH2F.h"

#include "plotting_functions.h"

void create_and_print_wf_graph_lines2 (std::string filename, std::vector<int> time, std::shared_ptr<std::vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2, std::vector<std::pair<int,double>> peaks) {

    TGraph *gWaveform = new TGraph();
    std::string newname = filename + "";

    TLine * line1 = new TLine(start1,level1,end1,level1);
    TLine * line2 = new TLine(start2,level2,end2,level2);

    TMarker *marker1 = new TMarker(peaks[0].first, peaks[0].second, 43); // Point (6, 36) with style 20
    TMarker *marker2 = new TMarker(peaks[1].first, peaks[1].second, 43); // Point (6, 36) with style 20

    for (int k = 0; k < time.size(); k++){

        gWaveform -> SetPoint ( k, time[k], (*ampl)[k]);
    }
    print_graph_lines2(gWaveform, newname, "Sample [#]", "ADC counts [#]", 0, 4000, line1, line2, marker1, marker2);
    
}

void print_graph_lines2 (TGraph *graph, std::string title, std::string x_axis, std::string y_axis, double yMin, double yMax, TLine *l1, TLine *l2, TMarker *p1, TMarker *p2){

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

void print_BAT_CAM_match( TH2F *h_track, std::vector<std::pair<double,double>> points, std::string mode, std::string title) {
    
    TCanvas* cmatch = new TCanvas(title.c_str(),title.c_str(), 700, 700);
    cmatch->cd();

    h_track->SetTitle(title.c_str());
    h_track->Draw("COLZ");

    for (const auto& point : points) {

        TMarker* marker = new TMarker(point.first, point.second, 0);

        if (mode == "bat") {
            marker->SetMarkerStyle(29);
            marker->SetMarkerColor(kRed);
        } 
        else if (mode == "cam") {
            marker->SetMarkerStyle(34);
            marker->SetMarkerColor(kBlack);
        }

        marker->SetMarkerSize(2);
        marker->Draw("same");
    }

    cmatch->DrawClone();
    cmatch->Write("BAT-PMT Match",TObject::kWriteDelete);
    delete cmatch;
}

void printTrackProfiles ( TH1D *h1, TH1D *h2, std::string title) {

    TCanvas* c_profiles = new TCanvas("c_profiles","c_profiles",1000,500); 
    c_profiles->Divide(1,2);

    c_profiles->cd(1);
    h1->SetTitle("Transversal Profile");
    h1->Draw();

    c_profiles->cd(2); 
    h2->SetTitle("Longitudinal Profile"); 
    h2->Draw();
    
    c_profiles->DrawClone();
    c_profiles->Write("Track profiles",TObject::kWriteDelete);
    delete c_profiles;

}


////////////////////////////// Older plotting functions ////////////////////////////////
/*
void print_histogram(TH1F *histo, std::string title, std::string x_axis, std::string y_axis){

    TCanvas *c = new TCanvas("","", 800, 500);
    c->cd();
    histo->SetTitle(title.c_str());
    histo->SetLineColor(kAzure-5);
    histo->SetLineWidth(1);
    histo->SetFillStyle(3003);
    histo->SetFillColor(kAzure+5);
    histo->GetYaxis()->SetTitle(y_axis.c_str());
    histo->GetYaxis()->SetTitleOffset(1.0);
    histo->GetYaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetTitleOffset(1.0);
    histo->GetXaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetTitle(x_axis.c_str());
    histo->Draw("hist");
    // histo->Write();
}

void print_graph_simple(TGraph *graph, std::string title, std::string x_axis, std::string y_axis){ //}, double yMin, double yMax){

    TCanvas *c = new TCanvas("","", 800, 500);
    c->cd();
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(x_axis.c_str());
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetXaxis()->SetTitleOffset(1);
    graph->GetYaxis()->SetTitle(y_axis.c_str());
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleOffset(1);
    graph->SetLineColor(kAzure-5);
    graph->SetMarkerColor(kAzure-5);
    // graph->SetMarkerStyle(2);
    // graph->GetYaxis()->SetRangeUser(0,2304);
    // graph->GetXaxis()->SetLimits(0,2304);
    graph->Draw("al");
} 

void create_and_print_wf_graph_simple (std::string filename, std::vector<int> time, shared_ptr<std::vector<double>> ampl, std::string tag) {

    TGraph *gWaveform = new TGraph();
    std::string newname = filename + "_" + tag;

    for (int k = 0; k < time.size(); k++){

        gWaveform -> SetPoint ( k, time[k], (*ampl)[k]);
    }
    print_graph_simple(gWaveform, newname, "t [ms]", "Amplitude [mV]");
    gWaveform->SetName(newname.c_str());
    // gWaveform->Write(newname.c_str(),TObject::kWriteDelete);
}

void create_and_print_wf_graph_lines (std::string filename, std::vector<int> time, shared_ptr<std::vector<double>> ampl, double start, double end, double level) {

    TGraph *gWaveform = new TGraph();
    std::string newname = filename + "";

    TLine * line1 = new TLine(start,level,end,level);

    for (int k = 0; k < time.size(); k++){

        gWaveform -> SetPoint ( k, time[k], (*ampl)[k]);
    }
    print_graph_lines(gWaveform, newname, "Sample [#]", "ADC counts [#]", 0, 4000, line1);
    
    // I need to save the canvas to also print the lines
    // gWaveform->SetName(newname.c_str());
    // gWaveform->Write(newname.c_str(),TObject::kWriteDelete);
}

void print_graph_lines (TGraph *graph, std::string title, std::string x_axis, std::string y_axis, double yMin, double yMax, TLine *l1){

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
    graph->SetLineWidth(2);
    graph->SetMarkerColor(kAzure-5);
    // graph->GetYaxis()->SetRangeUser(yMin,yMax);
    graph->Draw("apl");

    l1->SetLineColor(kRed-4);
    l1->SetLineWidth(3);
    l1->SetLineStyle(1);
    l1->Draw("same");

    c->SetName(title.c_str());
    c->Write(title.c_str(),TObject::kWriteDelete);
} 
*/