#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <experimental/filesystem>
#include <string>
#include <numeric>
// #include "Analyzer.h"
#include "Track_analyzer/Analyzer.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TPolyLine3D.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TGraphMultiErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TLine.h"
#include "TROOT.h"
#include "TKey.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TTree.h"
#include "TBrowser.h"

using namespace std;
namespace fs = std::experimental::filesystem;

vector <string> trim_name( string full_name, char delimiter); 
void print_histogram(TH1F *histo, string title, string x_axis, string y_axis);
void print_graph_simple(TGraph *graph, string title, string x_axis, string y_axis);
void create_and_print_wf_graph_simple (string filename, vector<int> time, shared_ptr<vector<double>> ampl, string tag);
void create_and_print_wf_graph_lines (string filename, vector<int> time, shared_ptr<vector<double>> ampl, double start, double end, double level);
void create_and_print_wf_graph_lines2 (string filename, vector<int> time, shared_ptr<vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2, double p1t, double p1a, double p2t, double p2a);
void print_graph_lines (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax, TLine *l1);
void print_graph_lines2 (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax, TLine *l1, TLine *l2, TMarker *p1, TMarker *p2);
void build_3D_vector (double x0, double x1, double y0, double y1, double z0, double z1,
 double l_xy, double a_xy, double l_z, int d_z, double p_z, double a_z, double length,
 int ru, int pi, int tr);

void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E);
void movingAverageFilter(std::shared_ptr<std::vector<double>>& input, int windowSize);
void create_bat_input(vector<vector<double>> slice_ints);
void run_then_read_bat(bool verbose, vector<double>& x ,vector<double>& y );
tuple <double,double> coord_change_pmt(double x, double y);

// struct Alpha3D {

//     int run;
//     int pic;
//     int trg;

//     double angle_XY; 
//     double trv_XY;
//     int quad;

//     double IP_X_cm;
//     double IP_Y_cm;
// };

struct AlphaTrackCAM {

    int run;
    int pic;

    double angle_XY; 
    double trv_XY;
    int quad;

    double IP_X_cm;
    double IP_Y_cm;

    bool its_alpha;
};

struct AlphaTrackPMT {

    int run;
    int pic;
    int trg;

    int dir; // -1 = towards GEM ; 1 = towards cathode; 0 = ambiguous
    double prob; // -1 = towards GEM ; 1 = towards cathode; 0 = ambiguous
    double trv_Z;
    int quad;

    bool its_alpha;

};

struct bayesd_data {

    double run, picture, trigger, peakIndex, L, L_std, x, x_std, y, y_std;
};
void read_file_fitted_results (string file, vector<bayesd_data> &data);


int main(int argc, char**argv) {
    std::cout << std::boolalpha; // Enable boolalpha
    TApplication *myapp=new TApplication("myapp",0,0); cout << "\n" << endl;

    string filename_cam = argv[ 1 ];

    string filename_pmt = argv[ 2 ]; 

    string outputfile = argv[ 3 ]; 
    TFile* file_root = new TFile(outputfile.c_str(),"recreate");

    // ofstream outFile("output_for_bat.txt");

    int debug_event = atof(argv[ 4 ]);
    // int debug_trigger = atof(argv[ 5 ]);


    /* ****************************************  Definition of variables and graphs  ******************************************************************  */

    // TH2F *hSpace         = new TH2F("","",2304,0,2305,2304,0,2305);

    vector<int> BeginScPix;
    vector<int> EndScPix;
    int nSc_red=0;

    vector<int> time_fast_wf(1024); // Create a vector with 1025 elements
    vector<int> time_slow_wf(4000); // Create a vector with 1025 elements

    iota(time_fast_wf.begin(), time_fast_wf.end(), 0); // Fill the vector with values starting from 0
    iota(time_slow_wf.begin(), time_slow_wf.end(), 0); // Fill the vector with values starting from 0

    int direction; // -1 = towards GEM ; 1 = towards cathode; 0 = ambiguous
    vector<double> total_integral;
    vector<int> skew_sign;
    vector<double> skew_ratio;          // ** check how to calculate skewness. Useful to then connect with Bayes fit.
    double avg_skew;

    double delta_z;
    vector<double> travelled_Z;
    double avg_travel_z;

    // Variables for ToT
    int TOT20_begin = 0, TOT20_end = 0;
    double TOT20_div;
    vector<double> TOT30;
    vector<double> TOT20;

    // vector<vector<double>> wfs(4);

    bool pmt_PID1 = false, pmt_PID2 = false;
    bool pmt_PID_total;
    int count_id1, count_id2;

    // electron drift velocity in LIME, He:CF4 60:40, 800 V/cm 
    const double drift_vel = 5.471; // cm/µs

    double first_half_int = 0, second_half_int = 0;
    double TOT20_half_way;
    int quadrant_pmt;


    vector<vector<double>> peak_int;

    vector<double> x_bat, y_bat;
    TH2F* h_fullsize = new TH2F(); 


    /* ************* CAM variables ********** */

    //impact point and directionality
    // For alphas (David)
    Int_t NPIP=2000;
    Float_t wFac=3.5;
  
    // For ER (original from Samuele)
    // Int_t NPIP=80;
    // Float_t wFac=2.5;

    double x_impact, y_impact;
    double xbar, ybar;
    int quadrant_cam;
    double angle_cam;
    const double granularity = 0.0155; // cm/pixel 

    bool cam_PID = false;

    /* ************* Combined information ********** */

    double x0, y0, z0;
    double x1, y1, z1;
    double theta_angle;

    /* ****************************************  Opening root recoed file -- CAMERA ******************************************************************  */

    TFile *reco_data_cam = TFile::Open(filename_cam.c_str());
    vector <string> name_trim1 = trim_name(filename_cam, '/');
    cout << "CAM Reco data file openend: " << name_trim1.back() << endl;

    TTree *tree_cam = (TTree*)reco_data_cam->Get("Events");

    int cam_run;        tree_cam->SetBranchAddress("run",   &cam_run);
    int cam_event;      tree_cam->SetBranchAddress("event", &cam_event);

    vector<float> sc_integral;    sc_integral.reserve(150000);    tree_cam->SetBranchAddress("sc_integral",   sc_integral.data());
    vector<float> sc_rms;         sc_rms.reserve(150000);         tree_cam->SetBranchAddress("sc_rms",        sc_rms.data());
    vector<float> sc_tgausssigma; sc_tgausssigma.reserve(150000); tree_cam->SetBranchAddress("sc_tgausssigma",sc_tgausssigma.data());
    vector<float> sc_length;      sc_length.reserve(150000);      tree_cam->SetBranchAddress("sc_length",     sc_length.data());
    vector<float> sc_width;       sc_width.reserve(150000);       tree_cam->SetBranchAddress("sc_width",      sc_width.data());
    vector<float> sc_xmean;       sc_xmean.reserve(150000);       tree_cam->SetBranchAddress("sc_xmean",      sc_xmean.data());  
    vector<float> sc_ymean;       sc_ymean.reserve(150000);       tree_cam->SetBranchAddress("sc_ymean",      sc_ymean.data());  
    vector<float> sc_nhits;       sc_nhits.reserve(150000);       tree_cam->SetBranchAddress("sc_nhits",      sc_nhits.data());  

    UInt_t nSc;                     tree_cam->SetBranchAddress("nSc",   &nSc);
    UInt_t Nredpix=0;               tree_cam->SetBranchAddress("nRedpix",&Nredpix);

    // vector<UInt_t> ScNpixels;       ScNpixels.reserve(nscmax);

    vector<float> sc_redpixID;      sc_redpixID.reserve(150000);    tree_cam->SetBranchAddress("sc_redpixIdx",sc_redpixID.data());
    vector<int> XPix;               XPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_iy",       XPix.data());  
    vector<int> YPix;               YPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_ix",       YPix.data());
    vector<float> ZPix;             ZPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_iz",       ZPix.data());  

    /* ****************************************  Opening root recoed file -- PMT ******************************************************************  */

    TFile *reco_data_pmt = TFile::Open(filename_pmt.c_str());
    vector <string> name_trim2 = trim_name(filename_pmt, '/');
    cout << "PMT Reco data file openend: " << name_trim2.back() << "\n" << endl;

    TTree *tree_pmt = (TTree*)reco_data_pmt->Get("PMT_Events");

    int pmt_run;        tree_pmt->SetBranchAddress("pmt_wf_run",        &pmt_run);
    int pmt_event;      tree_pmt->SetBranchAddress("pmt_wf_event",      &pmt_event);
    int pmt_trigger;    tree_pmt->SetBranchAddress("pmt_wf_trigger",    &pmt_trigger);
    int pmt_channel;    tree_pmt->SetBranchAddress("pmt_wf_channel",    &pmt_channel);
    int pmt_sampling;   tree_pmt->SetBranchAddress("pmt_wf_sampling",   &pmt_sampling);

    Float_t pmt_baseline;   tree_pmt->SetBranchAddress("pmt_wf_baseline",   &pmt_baseline);
    Float_t pmt_RMS;        tree_pmt->SetBranchAddress("pmt_wf_RMS",        &pmt_RMS);

    // Use this method to properly read the array since it's not a real vector and then otherwise reads wronmgly the addresses
    // The max length is 4000 (for slow). For reading the fast, I read just the first 1024 samples (the rest is "past memory")
    const int maxArraySize_fast = 1024;                 
    const int maxArraySize_slow = 4000;                  
    Float_t pmt_wf_ID[maxArraySize_slow];           tree_pmt->SetBranchAddress("pmt_wf_ID_full",pmt_wf_ID);
    Float_t pmt_waveforms[maxArraySize_slow];       tree_pmt->SetBranchAddress("pmt_fullWaveform_Y",pmt_waveforms);

    file_root->cd();

    /* **********************************************  Analysis  CAMERA  ********************************************** */

    vector<AlphaTrackCAM> CAM_alphas;
    Int_t nentries = (Int_t)tree_cam->GetEntries();
    for (Int_t cam_entry = 0; cam_entry < nentries; cam_entry++) {
        
        tree_cam->GetEntry(cam_entry);
        // cout << "Cam run: " << cam_run << "; event(pic): " << cam_event << "; nSc: " << nSc << endl;

        if ( cam_entry == debug_event ) {          // Choose specific event, for testing and debugging.

            file_root->mkdir(Form("Run_%i_ev_%i", cam_run, cam_event));
            file_root->cd(Form("Run_%i_ev_%i", cam_run, cam_event));

            sc_redpixID.clear();
            ScIndicesElem(nSc, Nredpix, sc_redpixID.data(), nSc_red, BeginScPix, EndScPix);

            for ( int sc_i = 0; sc_i < nSc_red; sc_i++ ) {

                cam_PID = false;
                cout << "\n\t*Cam run: " << cam_run << "; event: " << cam_event << "; cluster ID: " << sc_i << "\n" << endl;

                // if ( (sc_rms[nc] > 6) && (0.152 * sc_tgausssigma[nc] > 0.3) && (0.152 * sc_length[nc] < 80) && (sc_integral[nc] > 1000) && (sc_width[nc] / sc_length[nc] > 0.8)            //cam cut
                // && ( ( (sc_xmean[nc]-1152)*(sc_xmean[nc]-1152) + (sc_ymean[nc]-1152)*(sc_ymean[nc]-1152) ) <(800*800)) ) {
                if ( sc_integral[sc_i]/sc_nhits[sc_i] > 25 && sc_length[sc_i] > 100 && sc_width[sc_i] > 50 ) {   //Alpha cut fropm Giorgio
                    
                    cam_PID = true;
                    //----------- Build Analyser track  -----------//

                    Analyzer Track(Form("Track_run_%i_ev_%i", cam_run, cam_event),XPix.data(),YPix.data(),ZPix.data(),BeginScPix[sc_i],EndScPix[sc_i]);

                    Track.SetWScal(wFac);
                    Track.SetNPIP(NPIP);
                    Track.ApplyThr();
                    // Track.RemoveNoise(50);
                    Track.RemoveNoise(75);
                    // Track.RemoveNoise(100);
                    Track.ImpactPoint(Form("Track_event%i_run%i", cam_run, cam_event));
                    Track.ScaledTrack(Form("Track_event%i_run%i", cam_run, cam_event));
                    Track.Direction();
                    Track.ImprCorrectAngle();
                    Track.BuildLineDirection();

                    h_fullsize = Track.PlotandSavetoFileCOLZ_fullSize(Form("Track_event%i_run%i", cam_run, cam_event));


                    //----------- Get important track parameters  -----------//

                    angle_cam = Track.GetDir()/TMath::Pi()*180.;
                    xbar = Track.GetXbar();
                    ybar = Track.GetYbar();
                    x_impact = Track.GetXIP() * granularity;
                    y_impact = Track.GetYIP() * granularity;

                    // The y needs to be inverted for this calculation from the way root makes the plots.
                    // TO BE CROSS CHECKED WITH NEW DATA.
                    if      ( xbar < 2305./2. && ybar > 2305./2.) quadrant_cam = 1;
                    else if ( xbar > 2305./2. && ybar > 2305./2.) quadrant_cam = 2;
                    else if ( xbar > 2305./2. && ybar < 2305./2.) quadrant_cam = 3;
                    else if ( xbar < 2305./2. && ybar < 2305./2.) quadrant_cam = 4; 

                    //----------- Verbose information  -----------//

                    cout << "--> The particle in this cluster was identified as an alpha: " << cam_PID << endl;
                    cout << "\nTrack information: \n" << endl; 
                    cout << "--> Position barycenter: " << "x: " << xbar << "; y: " << ybar << endl;
                    cout << "--> Quadrant: " << quadrant_cam << endl;
                    cout << "--> Angle: " << angle_cam << " degrees." << endl;
                    cout << "--> Length (cm): " << sc_length[sc_i] * granularity << endl;

                    // Manually build the track in a TH2F. Not needed anymore.
                    /*for ( int pixel = BeginScPix[sc_i]; pixel < EndScPix[sc_i]; pixel++ ){
                        if ( ZPix[pixel] > 0 ) {
                            
                            // hSpace->SetBinContent(hSpace->GetXaxis()->FindBin(XPix[pixel]), hSpace->GetXaxis()->FindBin(YPix[pixel]), ZPix[pixel]);
                            // NOTA BENE: This correction makes sense for this recoed data. This should always be confirmed with a "debug" picture with the new version of the reco.
                            
                            // hSpace->SetBinContent(hSpace->GetXaxis()->FindBin(2305-XPix[pixel]), hSpace->GetXaxis()->FindBin(2305-YPix[pixel]), ZPix[pixel]);
                            hSpace->Fill(2305-XPix[pixel], 2305-YPix[pixel], ZPix[pixel]);      
                        }
                    }
                    TCanvas* c_track_cam = new TCanvas("c_track_cam","c_track_cam",3500,1500);
                    c_track_cam->cd();      
                    hSpace->DrawCopy("COLZ");
                    c_track_cam->Write("Picture",TObject::kWriteDelete);
                    */

                    //----------- Analyser saving options  -----------//

                    // Track.SavetoFile("test1");               // Saves TH2D in root file
                    // Track.SavetoFileDavid("test5");          // Saves TH2D, colz, in root file
                    // Track.SavePic("test2.png");              // saves png of track in folder with a bunch of variables.
                    // Track.SavePicDir("test3.png");           // Saves 3 plots with the real direction calculation with Samuele's code
                    // Track.SaveRootFile("test4.root");        // Saves a couple plots in a different root file.
                    Track.PlotandSavetoFileCOLZ(Form("Track_ev%i_run%i", cam_run, cam_event));            // Saves TH2D, colz, in root file
                    Track.PlotandSavetoFileDirectionalFull("X-Y Analyser");        // Saves the directionality plots all together, in root file

                    // Track profiles
                    TCanvas* c_profile = new TCanvas("c_profile","c_profile",1000,500); 
                    c_profile->Divide(1,2);

                    c_profile->cd(1);
                    TH1D* TrackProfileTrans = Track.FillProfile(false);
                    TrackProfileTrans->SetTitle("Transversal Profile");
                    TrackProfileTrans->Draw();

                    c_profile->cd(2); 
                    TH1D* TrackProfileLongi = Track.FillProfile(true);
                    TrackProfileLongi->SetTitle("Longitudinal Profile"); 
                    TrackProfileLongi->Draw();
                    c_profile->Write("Track profiles",TObject::kWriteDelete);
                    c_profile->DrawClone();
                    delete c_profile;

                    //----------- Collect all the relevant info for posterior analysis  -----------//

                    CAM_alphas.push_back({
                        .run = cam_run,
                        .pic  = cam_event,

                        .angle_XY = angle_cam,
                        .trv_XY = sc_length[sc_i] * granularity,

                        .quad = quadrant_cam,

                        .IP_X_cm = x_impact,
                        .IP_Y_cm = y_impact,

                        .its_alpha = cam_PID
                    });
                } else cout << "--> The particle in this cluster was identified as an alpha: " << cam_PID << endl;
            }
        }
        sc_redpixID.resize(nSc); // Resize to avoid problem with vector sizes
    }
    reco_data_cam->Close();
    

    /* *********************************************  Analysis PMT  *******************************************  */

    vector<AlphaTrackPMT> PMT_alphas;
    Int_t nentries_pmt = (Int_t)tree_pmt->GetEntries();
    for (Int_t pmt_wf = 0; pmt_wf < nentries_pmt; pmt_wf++) {  // Each entry is one waveform
        
        tree_pmt->GetEntry(pmt_wf);

        // cout << "PMT run: " << pmt_run << "; event: " << pmt_event << "; trigger: " << pmt_trigger << "; channel: " << pmt_channel << "; sampling: " << pmt_sampling << endl;

        shared_ptr<vector<double>> fast_waveform = make_shared<vector<double>>();
        shared_ptr<vector<double>> slow_waveform = make_shared<vector<double>>();

        // FOR PURE PMT ANALYSIS, A CUT TO SELECT ALPHAS WOULD BE NEEDED
        // TO DO MIXED ANALYSIS, I CAN AND SHOULD CUT ON THE CAMERA (TO ALLOW ASSOCIATION)

        // if ( pmt_event == debug_event && pmt_trigger == debug_trigger && pmt_sampling == 1024) {         // Choose specific event, for testing and debugging.
        // if ( pmt_event == debug_event && pmt_sampling == 1024) {         // Choose specific event, for testing and debugging.
        if ( pmt_event == debug_event && pmt_trigger == 1 && pmt_sampling == 1024) {         // Choose specific event, for testing and debugging.

            pmt_PID1 = false, pmt_PID2 = false, pmt_PID_total = false;

            if ( pmt_channel == 1) {
                cout << "\n\t*PMT run: " << pmt_run << "; event: " << pmt_event << "; trigger: " << pmt_trigger << "; sampling: " << pmt_sampling << endl;
            }

            //-----------  Put branch information into vector for further analysis  -----------//

            for (int j = 0; j < maxArraySize_fast; ++j) fast_waveform->push_back(pmt_waveforms[j]);

            //-----------  Moving average filter applied to smooth high frequency noisy  -----------//
            // Needed to properly get start and end. No problem in exaggerating as I'm not looking for any special features.

            movingAverageFilter(fast_waveform, 20);         

            //-----------  Getting max value of waveform to get ToT  -----------//

            auto max_it = max_element(fast_waveform->begin(), fast_waveform->end());        // Find the maximum value in the vector
            double max_value = *max_it;                                                     // Dereference the iterator to get the maximum value
            int time_max_value = time_fast_wf[(distance(fast_waveform->begin(), max_it))];

            TOT20_div = 0.20;                                                                 // Percentage of maximum used to get ToT level
            TOT20_begin = 0, TOT20_end = 0;
            for (int p = 0; p < fast_waveform->size(); ++p) {
                
                if ( ( (*fast_waveform)[p] > (max_value * TOT20_div) ) && TOT20_begin == 0)                                                         TOT20_begin   = time_fast_wf[p];
                if ( ( (*fast_waveform)[p] < (max_value * TOT20_div) ) && TOT20_begin !=0 && TOT20_end == 0 && ( (time_fast_wf[p] - TOT20_begin) > 10) )        TOT20_end     = time_fast_wf[p];

            }
            if (TOT20_end == 0) TOT20_end = 1024;
            // create_and_print_wf_graph_lines(Form("WF_run_%i_evt_%i_trg_%i_ch_%i",pmt_run,pmt_event,pmt_trigger,pmt_channel), time_fast_wf, fast_waveform, TOT20_begin, TOT20_end, max_value * TOT20_div);
            // cout << "**TOT20: " << (TOT20_end-TOT20_begin) << endl;
            TOT20.push_back( (TOT20_end-TOT20_begin) );

            double TOT30_div = 0.30;                                                                 // Percentage of maximum used to get ToT level
            double TOT30_begin = 0, TOT30_end = 0;
            for (int p = 0; p < fast_waveform->size(); ++p) {
                
                if ( ( (*fast_waveform)[p] > (max_value * TOT30_div) ) && TOT30_begin == 0)                     TOT30_begin   = time_fast_wf[p];
                if ( ( (*fast_waveform)[p] < (max_value * TOT30_div) ) && TOT30_begin != 0 && TOT30_end == 0 && ( (time_fast_wf[p] - TOT30_begin) > 10) )    TOT30_end     = time_fast_wf[p];
            }  
            if (TOT30_end == 0) TOT30_end = 1024;
            TOT30.push_back( (TOT30_end-TOT30_begin) );

            double ratio_pid = ( TOT20_end - TOT20_begin) / ( TOT30_end - TOT30_begin );

            // create_and_print_wf_graph_lines(Form("WF_run_%i_evt_%i_trg_%i_ch_%i",pmt_run,pmt_event,pmt_trigger,pmt_channel), time_fast_wf, fast_waveform, TOT30_begin, TOT30_end, max_value * TOT30_div);

            //-----------  Getting information for BAT  -----------//

            int n_slices = 5;
            double slice_width = (TOT20_end - TOT20_begin)/n_slices;
            double peak_int_c;
            double peak_val;

            peak_int.push_back(vector<double>{});
            for (int s = 0; s < n_slices; ++s) {

                peak_val = 0;
                peak_int_c = 0;
                for (int p = TOT20_begin+slice_width*(s); p < TOT20_begin+slice_width*(s+1); ++p) {
                
                    peak_val += (*fast_waveform)[p];
                    peak_int_c++;
                }
                peak_int[pmt_channel-1].push_back(peak_val/peak_int_c);
            }

            //-----------  Getting travelled Z  -----------//
            //    **SHOULD CORRECT FOR MINIMUM TILTNESS OFFSET
            
            delta_z = (TOT20_end - TOT20_begin) * (4/3) * (drift_vel/1000.); // in cm,         
            travelled_Z.push_back(delta_z);

            //-----------  Determination of direction (GEMs or cathode)  -----------//

            TOT20_half_way = (TOT20_begin + ((TOT20_end-TOT20_begin)/2.));
            first_half_int = 0, second_half_int = 0;
            for (int p_i = 0; p_i < fast_waveform->size(); ++p_i) {

                if ( p_i >= TOT20_begin && p_i < TOT20_half_way ) first_half_int  += (*fast_waveform)[p_i];
                if ( p_i > TOT20_half_way && p_i <= TOT20_end   ) second_half_int += (*fast_waveform)[p_i];
            }
            // cout << "first integral: " << first_half_int << "\nsecond integral: " << second_half_int << endl; 
            // cout << "ratio: " << (double)first_half_int/second_half_int << endl; 

            total_integral.push_back((double)(first_half_int + second_half_int));


            //----------- Calculation of skewness of the Bragg Peak  -----------//
            double peak1_amp = 0, peak1_time = 0;
            double peak2_amp = 0, peak2_time = 0;

            // if ( first_half_int > second_half_int ) {        // This doesn't always work because sometimes the max value is in the half with the *lowest* integral
            if ( time_max_value >= TOT20_begin && time_max_value < TOT20_half_way) {     // With this, I look for the "smaller peak" on the opposite half of the max value

                for (int p_s = TOT20_end - 1; p_s > TOT20_begin; --p_s) {

                    if ( (*fast_waveform)[p_s] < (*fast_waveform)[p_s+1] && (*fast_waveform)[p_s+1] > 0.4*max_value) { // this 0.4 cuts for peaks in the increasing slope

                        peak2_amp = (*fast_waveform)[p_s+1];
                        peak2_time = time_fast_wf[p_s+1];
                        peak1_amp = max_value;
                        peak1_time = time_max_value;
                        break;
                    } 
                }
            
                skew_sign.push_back(0);
                skew_ratio.push_back( (double)( ((-1.)*peak1_amp/peak2_amp) + 1. ) );
            }
            // else if ( second_half_int > first_half_int ) {
            else if ( time_max_value >= TOT20_half_way && time_max_value < TOT20_end ) {
                
                for (int p_s = TOT20_begin + 1; p_s < TOT20_end; ++p_s) {

                    if ( (*fast_waveform)[p_s] < (*fast_waveform)[p_s-1] && (*fast_waveform)[p_s-1] > 0.4*max_value ) {

                        peak1_amp = (*fast_waveform)[p_s-1];
                        peak1_time = time_fast_wf[p_s-1];
                        peak2_amp = max_value;
                        peak2_time = time_max_value;
                        break;
                    } 
                }

                skew_sign.push_back(1);
                skew_ratio.push_back( (double)( ((+1.)*peak2_amp/peak1_amp) - 1.) );
            }

            // cout << "Peak 1: " << peak1_time << " * " << peak1_amp << endl; 
            // cout << "Peak 2: " << peak2_time << " * " << peak2_amp << endl; 

            create_and_print_wf_graph_lines2(Form("WFtest_run_%i_evt_%i_trg_%i_ch_%i",pmt_run,pmt_event,pmt_trigger,pmt_channel), 
            time_fast_wf, fast_waveform,
            TOT20_begin, TOT20_end, max_value * TOT20_div, 
            TOT30_begin, TOT30_end, max_value * TOT30_div, 
            peak1_time, peak1_amp, peak2_time, peak2_amp);

            //----------- Calculation of all the final variables of interest  -----------//

            if ( pmt_channel == 4 ) {   // Perform actions only when I have the 4 waveforms

                bool verb_skew = false;
                if(verb_skew) {cout << "--> Skewness: " << endl;}                     
                if(verb_skew) {cout << "Ratios: ";}      
                if(verb_skew) {for (auto &val : skew_ratio) cout << val << " * "; cout<<endl;}

                auto max_skew_it = std::max_element(skew_ratio.begin(), skew_ratio.end(), [](double a, double b) { return abs(a) < abs(b);});

                double max_skew = 1;
                if (max_skew_it != skew_ratio.end()) {
    
                    max_skew = *max_skew_it; 
                    if(verb_skew) cout << "Abs max ratio: " << abs(max_skew) << endl;

                } else {

                    if(verb_skew) cout << "The vector is empty. This is not an alpha." << endl;
                    // skew_sign.clear();
                    // TOT30.clear();
                    // TOT20.clear();                
                    // skew_ratio.clear();                
                    // total_integral.clear();                
                    // travelled_Z.clear();
                    // continue;
                }

                if(verb_skew) {cout << "Normalized: ";}                     
                if(verb_skew) {for (auto &val : skew_ratio) cout << (double)val/abs(max_skew) << " # "; cout<<endl;}

                double avg_rat = 0;
                for (int elem = 0; elem < skew_ratio.size(); ++elem) {

                    if (elem != distance(skew_ratio.begin(), max_skew_it)) {

                        avg_rat += skew_ratio[elem]/abs(max_skew);  
                    } else if (elem == distance(skew_ratio.begin(), max_skew_it)) {

                        avg_rat += (max_skew/abs(max_skew)/2.);                            // Need to give some weight to the most clear wf. In this case, the weight is 1/2
                    }
                }

                if(verb_skew) {cout << "*Final probability: " << avg_rat*100. << " %" << endl;}

                cout << "\nPMT Track information: \n" << endl; 

                // DIRECTION -- To be further developed taking into account the real ratios, hopefully with the help of the bat fit.
                
                if      (avg_rat <= -0.5 )   direction = -1, cout << "--> Moving towards the GEMs with certainty: " <<  avg_rat*100. << " %" << endl;
                else if (avg_rat >= +0.5 )   direction = 1 , cout << "--> Moving towards the cathode with certainty: " << avg_rat*100. << " %" << endl;
                else if (avg_rat > -0.5 || avg_rat < +0.5 ) direction = 0 , cout << "--> Ambiguous. Certainty: " << avg_rat*100. << " %" << endl;
                

                /*  Old method 
                cout << "--> Skewness signs: "; for (auto i: skew_sign) cout << i << " ";

                avg_skew = accumulate(skew_sign.begin(), skew_sign.end(), 0.0) / skew_sign.size();
                cout << "  -> skewness sign average is: " << avg_skew << endl; // Related to how sure I am, I guess . Probably should use the real ratios (?)                  

                if      (avg_skew < 0.5 )   direction = -1, cout << "--> This track is moving towards the GEMs." << endl;
                else if (avg_skew >  0.5 )  direction = 1 , cout << "--> This track is moving towards the cathode." << endl;
                else if (avg_skew ==  0.5 ) direction = 0 , cout << "--> This track direction is ambiguous." << endl;
                // direction = round(avg_skew);
                */
                
                // TRAVELLED Z
                avg_travel_z = accumulate(travelled_Z.begin(), travelled_Z.end(), 0.0) / travelled_Z.size();
                cout << "--> The average travelled Z (cm) is: " << avg_travel_z << endl;

                // QUADRANT -- Could be improved by dividing the image in 16 pixels instead of 4.
                auto max_pmt_it = max_element(total_integral.begin(), total_integral.end()); 
                int index = distance(total_integral.begin(), max_pmt_it);
                quadrant_pmt = index + 1;
                cout << "--> The track is in the quadrant: " << quadrant_pmt << endl;

                // Particle ID

                count_id1 = 0;
                count_id2 = 0;
                for (int q = 0; q < TOT20.size(); ++q) {
                    
                    if ( (TOT20[q] / TOT30[q]) >= 1 && (TOT20[q] / TOT30[q]) <= 2 ) count_id1++;
                    if ( TOT20[q] >= 80 ) count_id2++;
                }

                if ( count_id1 == 4) pmt_PID1 = true;
                if ( count_id2 == 4) pmt_PID2 = true;
                
                if ( pmt_PID1 && pmt_PID2) {
                    pmt_PID_total = true;
                } else pmt_PID_total = false;

                cout << "--> The TOT20/TOT30 ratios were: " << (double)TOT20[0]/TOT30[0] << " " << (double)TOT20[1]/TOT30[1] << " " << (double)TOT20[2]/TOT30[2] << " " << (double)TOT20[3]/TOT30[3] << " " << endl;
                cout << "--> The TOT20 lengths were: " << TOT20[0] << " " << TOT20[1] << " " << TOT20[2] << " " << TOT20[3] << " " << endl;
                cout << "--> The particle in this trigger was identified as an alpha: " << pmt_PID_total << endl;

                ///////// Information for BAT
                for ( int pp = 0; pp < peak_int.size(); ++pp) {

                    cout << "Average values for slice integrals for pmt number: " << pp+1 << ": ";

                    for ( int qq = 0; qq < peak_int[pp].size(); ++qq ) {

                        cout << peak_int[pp][qq] << " ";
                        create_bat_input(peak_int);
                    }
                    cout << endl;
                }
                run_then_read_bat(true, x_bat, y_bat);


                // for (int i =0; i<x_bat.size();++i) {

                //     cout << x_bat[i]<< endl;
                // }

                TCanvas* cmatch = new TCanvas("cmatch", "cmatch", 700, 700);
                cmatch->cd();
                h_fullsize->Draw("");
                h_fullsize->SetName("BAT-PMT Match");
                for (int i =0; i<x_bat.size();++i) {

                    cout << x_bat[i] << " * ";
                    // Step 3: Add points using TMarker
                    TMarker* marker5 = new TMarker(x_bat[i], y_bat[i], 29);  // Cross marker
                    marker5->SetMarkerColor(kRed);
                    marker5->SetMarkerSize(2);
                    marker5->Draw("same");
                }
                cmatch->DrawClone();

                //----------- Collect all the relevant info for posterior analysis  -----------//

                PMT_alphas.push_back({
                    .run = pmt_run,
                    .pic  = pmt_event,
                    .trg = pmt_trigger,

                    .dir = direction,
                    .prob = (avg_rat*100.),
                    .trv_Z = avg_travel_z,

                    .quad = quadrant_pmt,

                    .its_alpha = pmt_PID_total
                });

                //----------- Clear vector used for the 4 waveforms  -----------//

                skew_sign.clear();
                TOT30.clear();
                TOT20.clear();                
                skew_ratio.clear();                
                total_integral.clear();                
                travelled_Z.clear();
                peak_int.clear();
                x_bat.clear();
                y_bat.clear();
            }
        }
    }
    reco_data_pmt->Close();

    /* **************************************  COMBINED ANALYSIS  ********************************************************************  */

    for (const auto& cam : CAM_alphas) {

        for (const auto& pmt : PMT_alphas) {
        
            if ( pmt.run == cam.run ) {

                if ( pmt.pic == cam.pic ) {
                
                    if (pmt.quad == cam.quad ) {        // Basic quadrant association

                        if ( pmt.its_alpha == true && cam.its_alpha == true) {

                            cout << "\n# Matched alpha in quadrant: " << pmt.quad << "; in trigger: " << pmt.trg << "; with Alpha-PID = " << pmt.its_alpha << endl;
                            
                            //----------- Creating 3D alpha track information  -----------//

                            x0 = cam.IP_X_cm;
                            y0 = cam.IP_Y_cm;               
                            z0 = 25.0;                          // Absolute Z fixed.

                            x1 = x0 + ( cam.trv_XY * cos(cam.angle_XY * TMath::Pi()/180.));
                            y1 = y0 + ( cam.trv_XY * sin(cam.angle_XY * TMath::Pi()/180.));
                            z1 = z0 + ( pmt.trv_Z  * pmt.dir);

                            double theta_angle = atan(pmt.trv_Z/cam.trv_XY) * 180. / TMath::Pi();

                            double length = TMath::Sqrt(pow(cam.trv_XY,2) + pow(pmt.trv_Z,2));

                            cout << "\n\t ** 3D Alpha track information: ** \n" << endl; 
                            cout << "--> Position, X: " << x0 << "; Y: " << y0  << endl;
                            cout << "--> Travelled XY: "    << cam.trv_XY   << endl;
                            cout << "--> Angle XY (#phi): "        << cam.angle_XY << endl;
                            cout << "--> Travelled Z: "     << pmt.trv_Z    << endl;
                            cout << "--> Direction in Z: "  << pmt.dir      << " at " << abs(pmt.prob) << " score" << endl;
                            cout << "--> Angle Z (#theta): "     << theta_angle  << endl;
                            cout << "--> 3D alpha length (cm): " << length  << endl;

                            build_3D_vector(x0,x1,y0,y1,z0,z1,cam.trv_XY, cam.angle_XY, pmt.trv_Z, pmt.dir, pmt.prob, theta_angle, length, pmt.run, pmt.pic, pmt.trg);
                        
                        } else continue;
                    } else continue;
                } else continue;
            } else continue;
        } //close for pmt
    } //close for cam


    file_root->Close();
    cout << "**Finished**" << endl;
    myapp->Run();

    return 0;
}


/*  ***************************************************************************************************************************  */
/*  ***************************************************************************************************************************  */
/*  ***************************************************************************************************************************  */

void create_bat_input(vector<vector<double>> slice_ints) {
    
    /*
    - **run**: The run number.
    - **event**: The event number.
    - **trigger**: The trigger number.
    - **peak index**: The index indicating the position of the peak in the waveform.
    - **L1**: The integral of the **PMT 1**.
    - **L2**: The integral of the **PMT 2**.
    - **L3**: The integral of the **PMT 3**.
    - **L4**: The integral of the **PMT 4**.

    --> The integral must be given in nC
    --> I[nC] = I[ADU] / 4096 * 4/3 * 1/50
    */
    
    double vtg_to_nC = (1./4096.) * (4./3.) * (1./50.);

    ofstream outFile("output_for_bat.txt");

    // Check if the file is open
    if (outFile.is_open()) {

        // Write values separated by tabs
        // for (int pmt_i = 0; pmt_i < slice_ints.size(); ++i) {

        for ( int slice_i = 0; slice_i < slice_ints[0].size(); ++slice_i ) {

            outFile << "20000" << "\t" << 1 << "\t" << 2 << "\t" << slice_i << "\t" << slice_ints[0][slice_i] *vtg_to_nC  << "\t" << slice_ints[1][slice_i] *vtg_to_nC << "\t" << slice_ints[2][slice_i] *vtg_to_nC << "\t" << slice_ints[3][slice_i] *vtg_to_nC << "\n";
        }

            // outFile << values[i];
            // if (i < values.size() - 1) {
            //     outFile << '\t'; // Add a tab after each value except the last one
            // }
        // }

        // Close the file
        outFile.close();
    }


}

void run_then_read_bat(bool verbose, vector<double> &x ,vector<double> &y ) {

    // Define the command to run the Bash script
    const char* command = "./BAT_PMTs/./runfit.out -i ./output_for_bat.txt -o bat_alpha_pos.txt -s 0 -e 50 -m association";

    // Run the command using the system function
    int result = system(command);

    // Check the result of the system call
    if (result == 0) {
        std::cout << "Script executed successfully." << std::endl;
    } else {
        std::cerr << "Error executing script." << std::endl;
    }

    double x_corr, y_corr;
    vector<double> x_bat, y_bat;

    vector<bayesd_data> fitted_data;

    string fitted_results = "bat_alpha_pos.txt";

    read_file_fitted_results(fitted_results, fitted_data);

    for ( int integral = 0; integral < fitted_data.size(); ++integral) {

        if (verbose) {
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
        }

        tie(x_corr,y_corr) = coord_change_pmt(fitted_data[integral].x, fitted_data[integral].y);
        x.push_back(x_corr);
        y.push_back(y_corr);
    }
}

void read_file_fitted_results (string file, vector<bayesd_data> &data){
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





// Function to calculate indices for clusters based on reduced pixels
// nSc: number of superclusters
// npix: total number of pixels
// sc_redpixID: array of indices for the first reduced pixel in each cluster
// nSc_red: number of processed superclusters (output)
// B: begin indices for reduced pixels (output)
// E: end indices for reduced pixels (output)
void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E){
    nSc_red=0;
    B.clear();
    E.clear();

    vector<float> sc_redpix_start;
    sc_redpix_start.push_back(0); // Initialize the start indices

    int parcount=0; // Counter for iterating through pixels

    // Loop through each supercluster to process its reduced pixels
    for(int i=0; i<nSc; i++){
        if(sc_redpixID[i]>0){
            sc_redpix_start.push_back(sc_redpixID[i]); // Add the start index of the reduced pixels
        }
    }

    nSc_red = sc_redpix_start.size(); // Number of processed superclusters

    sc_redpix_start.push_back(npix); // Ensure the last pixel is included

    // Calculate the begin and end indices for each cluster's reduced pixels
    for(int i=0;i<sc_redpix_start.size()-1;i++){
        B.push_back(sc_redpix_start[i]);
        E.push_back(sc_redpix_start[i+1]);
    }

    sc_redpix_start.clear(); // Clear the temporary storage
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

void print_graph_simple(TGraph *graph, string title, string x_axis, string y_axis){ //}, double yMin, double yMax){

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

void create_and_print_wf_graph_simple (string filename, vector<int> time, shared_ptr<vector<double>> ampl, string tag) {

    TGraph *gWaveform = new TGraph();
    string newname = filename + "_" + tag;

    for (int k = 0; k < time.size(); k++){

        gWaveform -> SetPoint ( k, time[k], (*ampl)[k]);
    }
    print_graph_simple(gWaveform, newname, "t [ms]", "Amplitude [mV]");
    gWaveform->SetName(newname.c_str());
    // gWaveform->Write(newname.c_str(),TObject::kWriteDelete);
}

void create_and_print_wf_graph_lines (string filename, vector<int> time, shared_ptr<vector<double>> ampl, double start, double end, double level) {

    TGraph *gWaveform = new TGraph();
    string newname = filename + "";

    TLine * line1 = new TLine(start,level,end,level);

    for (int k = 0; k < time.size(); k++){

        gWaveform -> SetPoint ( k, time[k], (*ampl)[k]);
    }
    print_graph_lines(gWaveform, newname, "Sample [#]", "ADC counts [#]", 0, 4000, line1);
    
    // I need to save the canvas to also print the lines
    // gWaveform->SetName(newname.c_str());
    // gWaveform->Write(newname.c_str(),TObject::kWriteDelete);
}

void create_and_print_wf_graph_lines2 (string filename, vector<int> time, shared_ptr<vector<double>> ampl, double start1, double end1, double level1, double start2, double end2, double level2, double p1t, double p1a, double p2t, double p2a) {

    TGraph *gWaveform = new TGraph();
    string newname = filename + "";

    TLine * line1 = new TLine(start1,level1,end1,level1);
    TLine * line2 = new TLine(start2,level2,end2,level2);

    TMarker *marker1 = new TMarker(p1t, p1a, 43); // Point (6, 36) with style 20
    TMarker *marker2 = new TMarker(p2t, p2a, 43); // Point (6, 36) with style 20

    for (int k = 0; k < time.size(); k++){

        gWaveform -> SetPoint ( k, time[k], (*ampl)[k]);
    }
    print_graph_lines2(gWaveform, newname, "Sample [#]", "ADC counts [#]", 0, 4000, line1, line2, marker1, marker2);
    
    // I need to save the canvas to also print the lines
    // gWaveform->SetName(newname.c_str());
    // gWaveform->Write(newname.c_str(),TObject::kWriteDelete);
}

void print_graph_lines (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax, TLine *l1){

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

void movingAverageFilter(std::shared_ptr<std::vector<double>>& input, int windowSize) {
    // Ensure the window size is valid
    if (windowSize <= 0) {
        throw std::invalid_argument("Window size must be positive");
    }

    // Get the size of the input vector
    int n = input->size();
    
    // Create a temporary vector to hold the filtered values
    std::vector<double> filtered(n);

    // Apply the moving average filter
    for (int i = 0; i < n; ++i) {
        // Determine the window boundaries
        int start = std::max(0, i - windowSize / 2);
        int end = std::min(n, i + windowSize / 2 + 1);
        
        // Calculate the sum of the values in the window
        double sum = std::accumulate(input->begin() + start, input->begin() + end, 0.0);
        
        // Calculate the average and store it in the temporary vector
        filtered[i] = sum / (end - start);
    }

    // Copy the filtered values back to the original vector
    *input = std::move(filtered);
}

void build_3D_vector (double x0, double x1, double y0, double y1, double z0, double z1,
 double l_xy, double a_xy, double l_z, int d_z, double p_z, double a_z, double length,
 int ru, int pi, int tr){

    ///////////////////   For plotting purposes, Y and Z are swapped.

    //-------  3D box with LIME's dimensions in cm  --------//

    TCanvas *c_3D = new TCanvas("c_3D", Form("Alpha 3D vector run_%i_pic_%i_trig_%i", ru, pi, tr), 700, 700); c_3D->cd();

    TH3F *axis = new TH3F("h3", Form("Alpha 3D vector run_%i_pic_%i_trig_%i", ru, pi, tr) , 1, 0, 36, 1, 0, 50, 1, 0, 36);       //155, 36  
    // TH3F *axis = new TH3F("h3", Form("Alpha 3D vector run_%i_pic_%i_trig_%i", ru, pi, tr) , 1, -1.5, 34.5, 1, 0, 50, 1, -1.5, 34.5);   //155, 36    
    // TH3F *axis = new TH3F("h3", Form("Alpha 3D vector run_%i_pic_%i_trig_%i", ru, pi, tr) , 1, -0.9, 33.9, 1, 0, 50, 1, -0.9, 33.9);   //151, 34.8    
    axis->SetStats(0);
    axis->GetXaxis()->SetTitle("X");
    axis->GetYaxis()->SetTitle("Z");
    axis->GetZaxis()->SetTitle("Y");
    axis->Draw(""); // Need to use "Copy" to visualize the plot and also save the canvas

    //--------  Main 3D vector line  --------//

    /*  // To create a monocolor simple line
    double x[2] = {x0, x1};
    double z[2] = {y0, y1};
    double y[2] = {z0, z1};

    TPolyLine3D *line = new TPolyLine3D(2, x, y, z);
    line->SetLineColor(kAzure-5);
    line->SetLineWidth(4);
    line->Draw("same");
    */

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

    //--------  Create legend with other important information  --------//

    TLegend* l = new TLegend(0.60, 0.6, 0.95, 0.9);
    l->SetHeader("3D Alpha information", "L");
    l->SetTextAlign(12); // Align text left-top (vertical center)
    // l->SetMargin(0.05);  // Reduce margin to minimum
    l->AddEntry((TObject*)0,Form("Travelled XY = %.2f cm",              l_xy), "p");
    l->AddEntry((TObject*)0,Form("Angle XY (#phi) = %.2f #circ",        a_xy), "p");
    l->AddEntry((TObject*)0,Form("Travelled Z = %.2f cm",               l_z), "p");
    l->AddEntry((TObject*)0,Form("Direction in Z = %i at %.1f score",         d_z, abs(p_z)), "p");
    l->AddEntry((TObject*)0,Form("Angle Z (#theta) = %.2f #circ",       a_z), "p");
    l->AddEntry((TObject*)0,Form("3D alpha length (cm) = %.2f",   length), "p");
    l->Draw("same");
    c_3D->DrawClone();
    c_3D->Write("3D vector",TObject::kWriteDelete);
    delete c_3D;
} 