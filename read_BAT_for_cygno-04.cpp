#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TKey.h"

using namespace std;

struct bayesd_data {

    double run, picture, trigger, peakIndex, L, L_std, x, x_std, y, y_std;
};

class point {
	private:
		double X, Y;
	public:
		point () {X=0.; Y=0.;}
        point(double x, double y) : X(x), Y(y) {}

		void setPoint(double a, double b) {
			X = a;
			Y = b;
		}
		double getX(void) {
			return X;
		}
		double getY(void) {
			return Y;
		} 
};

vector <string> trim_name( string full_name, char delimiter); 
void read_file_fited_results (string file, vector<bayesd_data> &data);
tuple <double,double> coord_change_pmt(double x, double y);
tuple <double,double> coord_change_camera(double x, double y);
void print_histogram(TH1F *histo, string title, string x_axis, string y_axis);
void print_graph(TGraph *graph, string title, string x_axis, string y_axis);

int main(int argc, char**argv) {
    TApplication *myapp=new TApplication("myapp",0,0); cout << "\n" << endl;

    vector< vector< double> > pmt_events;

    string fitted_results = argv[ 1 ];

    /* ****************************************  Definition of variables and graphs  ******************************************************************  */

    vector<point> selected_pmt_clusters;

    double pmt_x_new, pmt_y_new;
    point p_pmt;

    TH1F *hIntegral         = new TH1F("","",100,0,40000);
    TH1F *hAllDistances_X   = new TH1F("","",100,-2000,2000);
    TH1F *hAllDistances_Y   = new TH1F("","",100,-2000,2000);
    TH1F *hAllDistances_R   = new TH1F("","",100,-100,2000);

    TGraph *gPMT_cls   = new TGraph(); int pp = 0;
    TGraph *gPMT_cls_selected   = new TGraph(); int qq = 0;

    double x_corr, y_corr;

    vector<bayesd_data> fitted_data;
    read_file_fited_results(fitted_results, fitted_data);


    for ( int integral = 0; integral < fitted_data.size(); ++integral) {

        cout << "*Info from this file:" << endl;
        cout << "run: " << fitted_data[integral].run << " * ";
        cout << "picture: " << fitted_data[integral].picture << " * ";
        cout << "trigger: " << fitted_data[integral].trigger << " * ";
        cout << "peakIndex: " << fitted_data[integral].peakIndex << " * ";
        cout << "L: " << fitted_data[integral].L << " * ";
        cout << "L_std: " << fitted_data[integral].L_std << " * ";
        cout << "x: " << fitted_data[integral].x << " * ";
        cout << "x_std: " << fitted_data[integral].x_std << " * ";
        cout << "y: " << fitted_data[integral].y << " * ";
        cout << "y_std: " << fitted_data[integral].y_std << " * ";
        cout << endl; // To separate entries for different `integral` values

        // tie(x_corr,y_corr) = coord_change_pmt(fitted_data[integral].x, fitted_data[integral].y);
        // gPMT_cls->SetPoint(pp, x_corr, y_corr); pp++;

        gPMT_cls->SetPoint(pp, fitted_data[integral].x, fitted_data[integral].y); pp++;
    }
    print_graph(gPMT_cls, "Fitted points positions", "X", "Y");

    myapp->Run();
    return 0;
}


void read_file_fited_results (string file, vector<bayesd_data> &data){
    // run_number[0]  event[1]  trigger[2]  peakIndex[3]  L[4]  L_std[5] x[6]  x_std[7]  y[8]  y_std[9]

    string line;
    ifstream myfile;

    double run, event, trigger, peakIndex, L, L_std, x, x_std, y, y_std;

    myfile.open(file.c_str());
    vector <string> name_trim = trim_name(file, '/');
    cout << "\nFitted results file opened: " << name_trim.back() << endl;
    cout << "Writing event to vector..." << endl;

    while ( getline ( myfile, line ) ) {

        istringstream iss ( line );	//creates string consisting of a line
        string token;

        getline (iss, token, '\t'); run = (double) stod(token);         
        getline (iss, token, '\t'); event = (double) stod(token);       
        getline (iss, token, '\t'); trigger = (double) stod(token);     
        getline (iss, token, '\t'); peakIndex = (double) stod(token);   
        getline (iss, token, '\t'); L = (double) stod(token);           
        getline (iss, token, '\t'); L_std = (double) stod(token);       
        getline (iss, token, '\t'); x = (double) stod(token);           
        getline (iss, token, '\t'); x_std = (double) stod(token);       
        getline (iss, token, '\t'); y = (double) stod(token);           
        getline (iss, token, '\t'); y_std = (double) stod(token);       

        data.push_back( {run, event, trigger, peakIndex, L, L_std, x, x_std, y, y_std} );
    }
    myfile.close();
}



tuple <double,double> coord_change_pmt(double x, double y) {   //pmt cm to pixels

    double x_new_tmp, y_new_tmp;
    x_new_tmp = x * 1970./33. + 180;
    // y_new_tmp = 2304 - (170 + y*1970./33.);
    y_new_tmp = y * 1970./33. + 370;

    return make_tuple( x_new_tmp, y_new_tmp );
}

tuple <double,double> coord_change_camera(double x, double y) {   //camera to conventional coordiantes

    double x_new_tmp, y_new_tmp;
    x_new_tmp = y;
    // y_new_tmp = 2304 - x;
    y_new_tmp = x;

    return make_tuple( x_new_tmp, y_new_tmp );
}

vector <string> trim_name( string full_name, char delimiter) {

    string trimmed_name;
    vector <string> tokens;
    stringstream check1(full_name);
    string intermediate;
    while(getline(check1, intermediate, delimiter)) tokens.push_back(intermediate);

    return tokens;
}

void print_histogram(TH1F *histo, string title, string x_axis, string y_axis){

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

void print_graph(TGraph *graph, string title, string x_axis, string y_axis){ //}, double yMin, double yMax){

    TCanvas *c = new TCanvas("","", 500, 800);
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
    graph->SetMarkerStyle(20);
    // graph->GetYaxis()->SetRangeUser(0,33);
    // graph->GetXaxis()->SetLimits(0,33);
    // graph->GetYaxis()->SetRangeUser(0,2304);
    // graph->GetXaxis()->SetLimits(0,2304);
    graph->GetYaxis()->SetRangeUser(0,80);
    graph->GetXaxis()->SetLimits(0,50);
    graph->Draw("ap");
} 