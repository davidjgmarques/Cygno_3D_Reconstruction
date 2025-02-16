#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TLine.h"
#include "TAxis.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TLatex.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TRandom3.h"
#include "TPaveStats.h"

#include "plotting_functions.h"

void create_and_print_wf_graph_lines3 (std::string filename, std::vector<int> time, std::shared_ptr<std::vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2, 
    std::vector<std::pair<int,double>> peaks_crown,  std::vector<std::pair<int,double>> peaks_energy_dep) {

    TGraph *gWaveform = new TGraph();
    std::string newname = filename + "";
    
    TLine * line1 = new TLine(start1,level1,end1,level1);
    TLine * line2 = new TLine(start2,level2,end2,level2);

    std::vector<TMarker*> markers_crown;
    for (const auto& peak : peaks_crown) {
        TMarker* marker1 = new TMarker(peak.first, peak.second, 43);
        markers_crown.push_back(marker1);
    }

    std::vector<TMarker*> markers_energy_dep;
    for (const auto& peak : peaks_energy_dep) {
        TMarker* marker2 = new TMarker(peak.first, peak.second, 43);
        markers_energy_dep.push_back(marker2);
    }

    for (int k = 0; k < time.size(); k++){
        gWaveform -> SetPoint ( k, time[k], (*ampl)[k]);
    }

    print_graph_lines3(gWaveform, newname, "Sample [#]", "ADC counts [#]", 0, 4000, line1, line2, markers_crown, markers_energy_dep);
}

void print_graph_lines3 (TGraph *graph, std::string title, std::string x_axis, std::string y_axis, double yMin, double yMax, TLine *l1, TLine *l2, std::vector<TMarker*> markers_cr, std::vector<TMarker*> markers_ed){

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

    for (const auto & marker1 : markers_cr) {
        marker1->SetMarkerColor(kOrange+10); 
        marker1->SetMarkerStyle(23);
        marker1->SetMarkerSize(1.5);
        marker1->Draw("same");
    }

    for (const auto & marker2 : markers_ed) {
        marker2->SetMarkerColor(kOrange); 
        marker2->SetMarkerStyle(23);
        marker2->SetMarkerSize(1.5);
        marker2->Draw("same");
    }

    c->SetName(title.c_str());
    c->Write(title.c_str());
}

void printTrackProfilesAndFit ( TH1D *h1, TH1D *h2, std::string title, TFitResultPtr &fitResult, bool save) {

    TF1* func = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
    func->SetParameters(1, h1->GetMean(), h1->GetRMS(), 0); 
    func->SetParName(0, "Amplitude");
    func->SetParName(1, "Mean");
    func->SetParName(2, "Sigma");
    func->SetParName(3, "Constant");

    // fitResult = h1->Fit(func, "RNSQ");
    fitResult = h1->Fit(func, "RNSE");

    TCanvas* c_profiles = new TCanvas("c_profiles","c_profiles",1000,800); 
    c_profiles->cd();
    c_profiles->Divide(1,2);

    c_profiles->cd(1);
    h1->SetTitle("Transversal Profile");
    h1->Draw();
    func->Draw("same");

    TPaveStats *pt = new TPaveStats(0.67, 0.67, 0.97, 0.97, "brNDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12); 
    pt->SetBorderSize(1); 
    pt->SetShadowColor(0);
    pt->SetTextFont(42);

    pt->AddText(Form("Entries        = %.0f", h1->GetEntries()));
    pt->AddText(Form("#chi^{2} / ndf = %.2f / %d",      func->GetChisquare(),  func->GetNDF()));
    pt->AddText(Form("C              = %.2f +/- %.2f",  func->GetParameter(0), func->GetParError(0)));
    pt->AddText(Form("#mu            = %.2f +/- %.2f",  func->GetParameter(1), func->GetParError(1)));
    pt->AddText(Form("#sigma         = %.2f +/- %.2f",  func->GetParameter(2), func->GetParError(2)));
    pt->AddText(Form("Baseline       = %.2f +/- %.2f",  func->GetParameter(3), func->GetParError(3)));
    pt->SetOptStat(1110);
    pt->SetOptFit(1110);
    pt->Draw();

    c_profiles->cd(2); 
    h2->SetTitle("Longitudinal Profile"); 
    h2->Draw();
    
    c_profiles->DrawClone();
    if (save) c_profiles->Write(Form("Track_profiles_%s", title.c_str()));
    delete c_profiles;
}

/*
void doubleGaussianTest( TH1D *h1, TH1D *h2, std::string title, TFitResultPtr &fitResult, bool save) {

    // int binmax = TransProfile_b->GetMaximumBin();
    // --> float center_b = TransProfile_b->GetXaxis()->GetBinCenter(binmax);
    // --> float std_b = TransProfile_b->GetRMS();
    // --> TF1* TransProfGauss_b = new TF1("TransProfGauss_b","gaus(0) + gaus(3)",0,TransProfile_b->GetNbinsX()-1);
    // --> TransProfGauss_b->SetParameters(1000, center_b,std_b, 100, 0, 10);
    // --> TransProfGauss_b->SetParameters(1000, center_b,std_b, 50, 0, 10);
    // --> TransProfGauss_b->SetParLimits(1,center_b-0.5,center_b+0.5);
    // --> TransProfile_b->Fit("TransProfGauss_b","QR","",center_b-3*std_b, center_b+3*std_b);}

    // TF1* dbGau = new TF1("dbGau", "[0]+ gaus(1)+gaus(4)", 0, h1->GetNbinsX()-1);
    TF1* dbGau = new TF1("dbGau", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2) + [6]", 0, h1->GetNbinsX()-1);
    dbGau->SetParameters(1000, h1->GetMean(),h1->GetRMS(),
                        100, 0, 10, 0);
    dbGau->SetParLimits(1,h1->GetMean()-0.5,h1->GetMean()+0.5);
    dbGau->SetParName(0, "Amplitude 1");
    dbGau->SetParName(1, "Mean 1");
    dbGau->SetParName(2, "Sigma 1");
    dbGau->SetParName(3, "Amplitude 2");
    dbGau->SetParName(4, "Mean 2");
    dbGau->SetParName(5, "Sigma 2");
    dbGau->SetParName(6, "Baseline");

    // dbGau = h1->Fit(dbGau, "RNSE");
    h1->Fit(dbGau, "RNSE");

    TCanvas* c_dbGau = new TCanvas("c_dbGau","c_dbGau",1000,800); 
    c_dbGau->cd();

    h1->SetTitle("Transversal double Gaus Profile");
    h1->Draw();
    dbGau->Draw("same");

    TPaveStats *pt = new TPaveStats(0.67, 0.67, 0.97, 0.97, "brNDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12); // Align text to the left
    pt->SetBorderSize(1); // Set border size to 1
    pt->SetShadowColor(0); // Remove shadow
    pt->SetTextFont(42); // Use a plain font
    
    pt->AddText(Form("Entries             = %.0f", h1->GetEntries()));
    pt->AddText(Form("#chi^{2} / ndf      = %.2f / %d",      dbGau->GetChisquare(),  dbGau->GetNDF()));
    pt->AddText(Form("C1                  = %.2f +/- %.2f",  dbGau->GetParameter(0), dbGau->GetParError(0)));
    pt->AddText(Form("#mu 1               = %.2f +/- %.2f",  dbGau->GetParameter(1), dbGau->GetParError(1)));
    pt->AddText(Form("#sigma 1            = %.2f +/- %.2f",  dbGau->GetParameter(2), dbGau->GetParError(2)));
    pt->AddText(Form("C2                  = %.2f +/- %.2f",  dbGau->GetParameter(3), dbGau->GetParError(3)));
    pt->AddText(Form("#mu 2               = %.2f +/- %.2f",  dbGau->GetParameter(4), dbGau->GetParError(4)));
    pt->AddText(Form("#sigma 2            = %.2f +/- %.2f",  dbGau->GetParameter(5), dbGau->GetParError(5)));
    pt->AddText(Form("Baseline            = %.2f +/- %.2f",  dbGau->GetParameter(6), dbGau->GetParError(6)));
    pt->SetOptStat(1110);
    pt->SetOptFit(1110);
    pt->Draw();

    c_dbGau->DrawClone();
    if (save) c_dbGau->Write(Form("Track_doubleGaus_%s", title.c_str()));
    delete c_dbGau;
}
*/

/* 
void centralGaussianTest( TH1D *h1, TH1D *h2, std::string title, TFitResultPtr &fitResult, bool save) {

    int binmax = h1->GetMaximumBin();
    float center_b = h1->GetXaxis()->GetBinCenter(binmax);
    float std_b = h1->GetRMS();

    // Ensure the fit is performed within the specified range
    double fitRangeMin = center_b - 0.5 * std_b;
    double fitRangeMax = center_b + 0.5 * std_b;
    std::cout << "Fit range: [" << fitRangeMin << ", " << fitRangeMax << "]" << std::endl;
    
    // TF1* centralGau = new TF1("centralGau", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0, h1->GetNbinsX()-1);
    TF1* centralGau = new TF1("centralGau", "[0]*exp(-0.5*((x-[1])/[2])^2)");
    centralGau->SetParameters(1000, h1->GetMean(),h1->GetRMS());
    // dbGau->SetParLimits(1,h1->GetMean()-0.5,h1->GetMean()+0.5);
    centralGau->SetParName(0, "Amplitude 1");
    centralGau->SetParName(1, "Mean 1");
    centralGau->SetParName(2, "Sigma 1");

    fitResult = h1->Fit(centralGau, "RS", "", fitRangeMin, fitRangeMax);

    TCanvas* c_centralGau = new TCanvas("c_centralGau","c_centralGau",1000,800); 
    c_centralGau->cd();

    h1->SetTitle("Transversal double Gaus Profile");
    h1->Draw();
    centralGau->Draw("same");

    TPaveStats *pt = new TPaveStats(0.67, 0.67, 0.97, 0.97, "brNDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12); // Align text to the left
    pt->SetBorderSize(1); // Set border size to 1
    pt->SetShadowColor(0); // Remove shadow
    pt->SetTextFont(42); // Use a plain font
    
    pt->AddText(Form("Entries             = %.0f", h1->GetEntries()));
    pt->AddText(Form("#chi^{2} / ndf      = %.2f / %d",      centralGau->GetChisquare(),  centralGau->GetNDF()));
    pt->AddText(Form("C1                  = %.2f +/- %.2f",  centralGau->GetParameter(0), centralGau->GetParError(0)));
    pt->AddText(Form("#mu 1               = %.2f +/- %.2f",  centralGau->GetParameter(1), centralGau->GetParError(1)));
    pt->AddText(Form("#sigma 1            = %.2f +/- %.2f",  centralGau->GetParameter(2), centralGau->GetParError(2)));
    pt->SetOptStat(1110);
    pt->SetOptFit(1110);
    pt->Draw();

    c_centralGau->DrawClone();
    if (save) c_centralGau->Write(Form("Track_centralGaus_%s", title.c_str()));
    delete c_centralGau;
} 
*/

double getHistogramRMS(TH1D* histo) {
    if (!histo) {
        std::cerr << "Error: Null histogram pointer provided." << std::endl;
        return -1;
    }
    return histo->GetRMS();
}

/*
Function to retrieve the length (extension) of the track profile above a given percentage of the maximum value
Possibly this function alrady exists in the form of other variables that I didn't know about
*/
double getTrackProfileWidth(TH1D* histo, double percentage = 0.05, int consecutiveBins = 3, bool verb = false) {
    if (!histo) {
        std::cerr << "Error: Null histogram pointer provided." << std::endl;
        return -1;
    }

    double max = histo->GetMaximum();
    double threshold = percentage * max;
    double width_prof = 0;
    bool aboveThreshold = false;
    double start = 0;
    double end = 0;
    int count = 0;

    for (int i = 1; i <= histo->GetNbinsX(); i++) {
        double binContent = histo->GetBinContent(i);
        double binCenter = histo->GetBinCenter(i);

        if (binContent >= threshold) {
            count++;
            if (count >= consecutiveBins) {
                if (!aboveThreshold) {
                    start = histo->GetBinCenter(i - consecutiveBins + 1); // Set start to the first of the consecutive bins
                    aboveThreshold = true;
                }
                end = binCenter;
            }
        } else {
            count = 0;
            if (aboveThreshold) {
                break;
            }
        }
    }

    if (aboveThreshold) {
        width_prof = end - start;
    }

    if (verb) std::cout << "Alpha track width from profile: " << width_prof << std::endl;

    return width_prof;
}

void Points_BAT_CAM(std::vector<std::pair<double,double>> batPoints, std::vector<std::pair<double,double>> camPoints, std::string title) {
    
    TCanvas *cpoints = new TCanvas("cpoints","cpoints", 800, 800);
    cpoints->cd();

    // Create a dummy histogram to set the axis ranges and titles
    TH2F *dummyHist = new TH2F("dummyHist", title.c_str(), 1, 0, 2304, 1, 0, 2304);
    dummyHist->GetXaxis()->SetTitle("X pixels [#]");
    dummyHist->GetXaxis()->SetTitleSize(0.045);
    dummyHist->GetXaxis()->SetTitleOffset(1);
    dummyHist->GetYaxis()->SetTitle("Y pixels [#]");
    dummyHist->GetYaxis()->SetTitleSize(0.045);
    dummyHist->GetYaxis()->SetTitleOffset(1);

    gStyle->SetOptStat(0);  // Disable statistics box
    dummyHist->Draw();

    TLegend *legend = new TLegend(0.6, 0.7, 0.87, 0.87);

    for (size_t i = 0; i < batPoints.size(); i++) {
        const auto& point = batPoints[i];
        TMarker* marker = new TMarker(point.first, point.second, 0);
        marker->SetMarkerStyle(29);
        marker->SetMarkerColor(kRed);
        marker->SetMarkerSize(2);
        marker->Draw("SAME");
        if (i == 0) legend->AddEntry(marker, "PMT", "p");
    }

    for (size_t i = 0; i < camPoints.size(); i++) {
        const auto& point = camPoints[i];
        TMarker* marker = new TMarker(point.first, point.second, 0);
        marker->SetMarkerStyle(34);
        marker->SetMarkerColor(kBlack);
        marker->SetMarkerSize(2);
        marker->Draw("SAME");
        if (i == 0) legend->AddEntry(marker, "CAM", "p");
    }

    legend->Draw("SAME");

    cpoints->DrawClone();
    cpoints->Write(title.c_str());
    delete cpoints;
    delete dummyHist; // Clean up the dummy histogram

}

void addTracks(TCanvas *image, TH2F* histo, TH2F* track, int fminx, int fminy, std::string nometh2) {

    int nBinsX = track->GetNbinsX();
    int nBinsY = track->GetNbinsY();

    int offsetX = fminx;
    int offsetY = fminy;

    // Iterate over the bins of the smaller histogram
    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            double content = track->GetBinContent(i, j);

            int binX = offsetX + static_cast<int>(i);
            int binY = offsetY + static_cast<int>(j);

            if (binX <= histo->GetNbinsX() && binY <= histo->GetNbinsY()) {
                histo->SetBinContent(binX, binY, content);
            }
        }
    }

    image->cd();
    histo->Draw("SAME COLZ");
    image->DrawClone();
    // canv->Write();
    
    // delete canv;
    // return histo;
}


void build_3D_vector (double x0, double x1, double y0, double y1, double z0, double z1,
 double l_xy, double a_xy, double l_z, int d_z, double p_z, double a_z, double length,
 int ru, int pi, int tr, bool cloud, double sigma){

    ///////////////////   For plotting purposes, Y and Z are swapped.

    //-------  3D box with LIME's dimensions in cm  --------//

    TCanvas *c_3D = new TCanvas(Form("Alpha_3D_vector_run_%i_pic_%i_trig_%i", ru, pi, tr), Form("Alpha_3D_vector_run_%i_pic_%i_trig_%i", ru, pi, tr), 700, 700); c_3D->cd();
    TH3F *axis = new TH3F(Form("Alpha_3D_vector_run_%i_pic_%i_trig_%i", ru, pi, tr), Form("Alpha_3D_vector_run_%i_pic_%i_trig_%i", ru, pi, tr) , 72, 0, 36, 100, 0, 50, 72, 0, 36); //nb: some bins are needed just to be able to zoom      
    axis->SetStats(0);
    axis->GetXaxis()->SetTitle("X");
    axis->GetYaxis()->SetTitle("Z");
    axis->GetZaxis()->SetTitle("Y");
    axis->Draw(""); // Need to use "Copy" to visualize the plot and also save the canvas

    //--------  Main 3D vector line  --------//

    // Normalize the vector for the arrowhead calculation
    double ux = (x1-x0);
    double uy = (y1-y0);
    double uz = (z1-z0);

    // To create a color gradient in the track - THIS PART COULD / SHOULD BE ASSOCIATED TO THE LONGITUDINAL PROFILE (?)
    const int numSegments = 45;
    const int color_jumper = 5;

    // Create arrays for the segment points
    double x[numSegments + 1];
    double y[numSegments + 1];
    double z[numSegments + 1];

    // Calculate the segment points and colors
    for (int i = 0; i <= numSegments; ++i) {
        double t = (double)i / numSegments;
        x[i] = x0 + t * ux;
        y[i] = y0 + t * uy;
        z[i] = z0 + t * uz;
    }

    // Create and draw the segments
    auto colorIndex = TColor::GetPalette();
    for (int i = 0; i < numSegments; ++i) {
        double segmentX[2] = {x[i], x[i+1]};
        double segmentZ[2] = {y[i], y[i+1]};
        double segmentY[2] = {z[i], z[i+1]};
        
        TPolyLine3D *segment = new TPolyLine3D(2, segmentX, segmentY, segmentZ);
        
        // Set color based on gradient factor
        segment->SetLineColor(colorIndex.At(i*color_jumper));
        segment->SetLineWidth(3);
        segment->Draw();
    }

    //--------  Create arrow lines to make direction more clear  --------//

    // Length of the arrowhead lines
    // double arrowHeadLength = 0.025 * length;
    double arrowHeadLength = 0.1;

    double x2 = x1 - arrowHeadLength * (ux + uy);
    double y2 = y1 - arrowHeadLength * (uy - ux);
    double z2 = z1 - arrowHeadLength * uz;

    double x3 = x1 - arrowHeadLength * (ux - uy);
    double y3 = y1 - arrowHeadLength * (uy + ux);
    double z3 = z1 - arrowHeadLength * uz;

    // Calculate orthogonal arrowhead points
    double x4 = x1 - arrowHeadLength * (ux + uz);
    double y4 = y1 - arrowHeadLength * uy;
    double z4 = z1 - arrowHeadLength * (uz - ux);

    double x5 = x1 - arrowHeadLength * (ux - uz);
    double y5 = y1 - arrowHeadLength * uy;
    double z5 = z1 - arrowHeadLength * (uz + ux);

    // Create arrays for the arrowhead lines
    double arrowX1[2] = {x1, x2};
    double arrowZ1[2] = {y1, y2};
    double arrowY1[2] = {z1, z2};

    double arrowX2[2] = {x1, x3};
    double arrowZ2[2] = {y1, y3};
    double arrowY2[2] = {z1, z3};

    double arrowX3[2] = {x1, x4};
    double arrowZ3[2] = {y1, y4};
    double arrowY3[2] = {z1, z4};

    double arrowX4[2] = {x1, x5};
    double arrowZ4[2] = {y1, y5};
    double arrowY4[2] = {z1, z5};

    // Create the arrowhead lines
    TPolyLine3D *arrowLine1 = new TPolyLine3D(2, arrowX1, arrowY1, arrowZ1);
    arrowLine1->SetLineColor(colorIndex.At(numSegments*color_jumper));
    arrowLine1->SetLineWidth(3);
    arrowLine1->Draw("same");

    TPolyLine3D *arrowLine2 = new TPolyLine3D(2, arrowX2, arrowY2, arrowZ2);
    arrowLine2->SetLineColor(colorIndex.At(numSegments*color_jumper));
    arrowLine2->SetLineWidth(3);
    arrowLine2->Draw("same");

    TPolyLine3D *arrowLine3 = new TPolyLine3D(2, arrowX3, arrowY3, arrowZ3);
    arrowLine3->SetLineColor(colorIndex.At(numSegments*color_jumper));
    arrowLine3->SetLineWidth(3);
    arrowLine3->Draw("same");

    TPolyLine3D *arrowLine4 = new TPolyLine3D(2, arrowX4, arrowY4, arrowZ4);
    arrowLine4->SetLineColor(colorIndex.At(numSegments*color_jumper));
    arrowLine4->SetLineWidth(3);
    arrowLine4->Draw("same");

    //--------  Electron cloud   ---------------------------------------//
    /*

    // Generate Gaussian points around the main line
    const int numGaussianPoints = 10000;
    double gaus_sigma = sigma;  // Standard deviation for the Gaussian distribution

    TRandom3 rand;
    int nbins = 100;
    TH3F *densityHist = new TH3F(Form("electroncloud_run_%i_pic_%i_trig_%i", ru, pi, tr),Form("electroncloud_run_%i_pic_%i_trig_%i", ru, pi, tr), nbins, 0, 36, nbins, 0, 50, nbins, 0, 36);

    for (int i = 0; i <= numSegments; ++i) {
        for (int j = 0; j < numGaussianPoints / numSegments; ++j) {
            double dx = rand.Gaus(0, gaus_sigma);
            double dy = rand.Gaus(0, gaus_sigma);
            double dz = rand.Gaus(0, gaus_sigma);

            double px = x[i] + dx;
            double pz = y[i] + dy;
            double py = z[i] + dz;

            densityHist->Fill(px, py, pz);  // Fill the histogram with the point
        }
    }

    densityHist->GetXaxis()->SetTitle("X");
    densityHist->GetYaxis()->SetTitle("Z");
    densityHist->GetZaxis()->SetTitle("Y");
    densityHist->Draw("SAME BOX2Z");  // Draw the histogram with color representation
    */

    //--------  Electron cloud   ---------------------------------------//

    if (cloud) {

        colorIndex = TColor::GetPalette();
        TRandom3 rand;
        const int numGaussianPoints = 2000;  //typical
        // const int numGaussianPoints = 200;  //for the gif
        double gaus_sigma = sigma;

        for (int i = 0; i <= numSegments; ++i) {
            for (int j = 0; j < numGaussianPoints / numSegments; ++j) {
                double dx = rand.Gaus(0, gaus_sigma);
                double dy = rand.Gaus(0, gaus_sigma);
                double dz = rand.Gaus(0, gaus_sigma);

                double px = x[i] + dx;
                double pz = y[i] + dy;
                double py = z[i] + dz;

                TPolyMarker3D *point = new TPolyMarker3D(1);
                point->SetPoint(0, px, py, pz);

                point->SetMarkerColor(colorIndex.At(i * color_jumper));
                point->SetMarkerStyle(20);  
                point->SetMarkerSize(0.5);
                point->Draw("same");
            }
        }
    }
    //--------  Create legend with other important information  --------//

    TLegend* l = new TLegend(0.60, 0.6, 0.95, 0.9);
    l->SetHeader("3D Alpha information", "L");
    l->SetTextAlign(12); // Align text left-top (vertical center)
    // l->SetMargin(0.05);  // Reduce margin to minimum
    l->AddEntry((TObject*)0,Form("Travelled XY = %.2f cm",              l_xy), "p");
    l->AddEntry((TObject*)0,Form("Angle XY (#phi) = %.2f #circ",        a_xy), "p");
    l->AddEntry((TObject*)0,Form("Travelled Z = %.2f cm",               l_z), "p");
    l->AddEntry((TObject*)0,Form("Direction in Z = %i at %.2f score",   d_z, p_z), "p");
    l->AddEntry((TObject*)0,Form("Angle Z (#theta) = %.2f #circ",       a_z), "p");
    l->AddEntry((TObject*)0,Form("3D alpha length (cm) = %.2f",   length), "p");
    l->Draw("same");
    c_3D->DrawClone();
    c_3D->Write();
    delete c_3D;
} 