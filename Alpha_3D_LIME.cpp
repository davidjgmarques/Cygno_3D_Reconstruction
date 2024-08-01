#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <experimental/filesystem>
#include <string>
#include <numeric>
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

#include "../Track_analyzer/Analyzer.h"
#include "./plotting_functions.h"
#include "./waveform_analyser.h"
#include "./bat_functions.h"


using namespace std;
namespace fs = std::experimental::filesystem;

void build_3D_vector (double x0, double x1, double y0, double y1, double z0, double z1,
 double l_xy, double a_xy, double l_z, int d_z, double p_z, double a_z, double length,
 int ru, int pi, int tr);

void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E);
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

    int dir;        // -1 = towards GEM ; 1 = towards cathode; 0 = ambiguous
    double prob;    // -1 = towards GEM ; 1 = towards cathode; 0 = ambiguous
    double trv_Z;
    int quad;

    bool its_alpha;

};

int main(int argc, char**argv) {
    std::cout << std::boolalpha; // Enable boolalpha
    TApplication *myapp=new TApplication("myapp",0,0); cout << "\n" << endl;

    string filename_cam = argv[ 1 ];

    string filename_pmt = argv[ 2 ];

    string outputfile = argv[ 3 ]; 
    TFile* file_root = new TFile(outputfile.c_str(),"recreate");

    int debug_event = atof(argv[ 4 ]);
    // int debug_trigger = atof(argv[ 5 ]);


    /* ****************************************  Definition of variables and graphs  ******************************************************************  */

    vector<int> BeginScPix;
    vector<int> EndScPix;
    int nSc_red=0;

    vector<int> time_fast_wf(1024); // Create a vector with 1025 elements
    vector<int> time_slow_wf(4000); // Create a vector with 1025 elements

    iota(time_fast_wf.begin(), time_fast_wf.end(), 0); // Fill the vector with values starting from 0
    iota(time_slow_wf.begin(), time_slow_wf.end(), 0); // Fill the vector with values starting from 0

    int direction; // -1 = towards GEM ; 1 = towards cathode; 0 = ambiguous
    double dir_score;
    bool verbose_dir_score = false;

    double wf_integral;
    vector<double> integrals_wfs;

    vector<double> skew_ratio;          // ** check how to calculate skewness. Useful to then connect with Bayes fit.
    double avg_skew;

    vector<pair<int,double>> wf_peaks;

    double delta_z;
    vector<double> travelled_Z;
    double avg_travel_z;

    // Variables for ToT
    int TOT20_begin = 0, TOT20_end = 0;
    double TOT20_div;

    int TOT30_begin = 0, TOT30_end = 0;
    double TOT30_div;
    vector<double> TOT30;
    vector<double> TOT20;

    // vector<vector<double>> wfs(4);

    bool pmt_PID_total;

    // electron drift velocity in LIME, He:CF4 60:40, 800 V/cm 
    const double drift_vel = 5.471; // cm/µs
    
    double TOT20_half_way;
    int quadrant_pmt;

    vector<vector<double>> integrals_slices;

    vector<pair<double,double>> points_bat;
    vector<pair<double,double>> points_cam;

    TH2F* h_fullsize = new TH2F(); 

    // Number of slices to be used for the track matching (CAM <-> BAT)
    int matching_slices = 5;


    /* ************* CAM variables ********** */

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
    vector<int> XPix;               XPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_ix",       XPix.data());  
    vector<int> YPix;               YPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_iy",       YPix.data());
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

    // Float_t pmt_baseline;   tree_pmt->SetBranchAddress("pmt_wf_baseline",   &pmt_baseline);
    // Float_t pmt_RMS;        tree_pmt->SetBranchAddress("pmt_wf_RMS",        &pmt_RMS);

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
                cout << "\t*Cam run: " << cam_run << "; event: " << cam_event << "; cluster ID: " << sc_i << "\n" << endl;

                // if ( (sc_rms[nc] > 6) && (0.152 * sc_tgausssigma[nc] > 0.3) && (0.152 * sc_length[nc] < 80) && (sc_integral[nc] > 1000) && (sc_width[nc] / sc_length[nc] > 0.8)            //cam cut
                // && ( ( (sc_xmean[nc]-1152)*(sc_xmean[nc]-1152) + (sc_ymean[nc]-1152)*(sc_ymean[nc]-1152) ) <(800*800)) ) {
                if ( sc_integral[sc_i]/sc_nhits[sc_i] > 25 && sc_length[sc_i] > 100 && sc_width[sc_i] > 50 ) {   //Alpha cut fropm Giorgio
                    
                    cam_PID = true;

                    //----------- Build Analyser track  -----------//

                    Analyzer Track(Form("Track_run_%i_ev_%i", cam_run, cam_event),XPix.data(),YPix.data(),ZPix.data(),BeginScPix[sc_i],EndScPix[sc_i]);

                    Track.SetWScal(wFac);
                    Track.SetNPIP(NPIP);
                    Track.ApplyThr();
                    // Track.RemoveNoise(75);
                    Track.RemoveNoise(100);
                    Track.ImpactPoint(Form("Track_event%i_run%i", cam_run, cam_event));
                    Track.ScaledTrack(Form("Track_rebinned_event%i_run%i", cam_run, cam_event));
                    Track.Direction();
                    Track.ImprCorrectAngle();
                    Track.BuildLineDirection();

                    h_fullsize = Track.PlotandSavetoFileCOLZ_fullSize(Form("Track_event%i_run%i", cam_run, cam_event));
                    points_cam = Track.GetLinePoints(matching_slices,"edges");
                    print_BAT_CAM_match(h_fullsize, points_cam, "cam", "CAM Match");

                    //----------- Get important track parameters  -----------//

                    angle_cam = Track.GetDir()/TMath::Pi()*180.;
                    xbar = Track.GetXbar();
                    ybar = Track.GetYbar();
                    x_impact = Track.GetXIP() * granularity;
                    y_impact = Track.GetYIP() * granularity;

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

                    Track.PlotandSavetoFileCOLZ(Form("Track_ev%i_run%i", cam_run, cam_event));          // Saves TH2D, colz, in root file
                    Track.PlotandSavetoFileDirectionalFull("X-Y Analyser");                             // Saves the directionality plots all together, in root file

                    printTrackProfiles( Track.FillProfile(false), Track.FillProfile(true), "Track profiles");

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

        shared_ptr<vector<double>> fast_waveform = make_shared<vector<double>>();
        shared_ptr<vector<double>> slow_waveform = make_shared<vector<double>>();

        // if ( pmt_event == debug_event && pmt_trigger == debug_trigger && pmt_sampling == 1024) {         // Choose specific event, for testing and debugging.
        // if ( pmt_event == debug_event && pmt_sampling == 1024) {         // Choose specific event, for testing and debugging.
        if ( pmt_event == debug_event && pmt_trigger == 0 && pmt_sampling == 1024) {         // Choose specific event, for testing and debugging.

            
            if ( pmt_channel == 1) cout << "\n\t*PMT run: " << pmt_run << "; event: " << pmt_event << "; trigger: " << pmt_trigger << "; sampling: " << pmt_sampling << endl;

            //-----------  Put branch information into vector for further analysis  --------------------------------------------//


            for (int j = 0; j < maxArraySize_fast; ++j) fast_waveform->push_back(pmt_waveforms[j]);


            //-----------  Moving average filter applied to smooth high frequency noisy  ----------------------------------------//
            // Needed to properly get start and end. No problem in exaggerating as I'm not looking for any special features.


            movingAverageFilter(fast_waveform, 20);         


            //-----------  Calculate maximum amplitude and get its time coordinate  -----------------------------------------------//


            auto max_it = max_element(fast_waveform->begin(), fast_waveform->end());                // Find the maximum value in the vector
            double max_value_a = *max_it;                                                           // Dereference the iterator to get the maximum value
            int max_value_t = time_fast_wf[(distance(fast_waveform->begin(), max_it))];


            //-----------  Calculate TOTs 20 and 30 using the value of waveform  ------------------------------------------------//


            TOT20_div = 0.20;                                                                   // Percentage of maximum used to get ToT level
            TOT20_begin = 0, TOT20_end = 0;

            TOT30_div = 0.30;                                                                   // Percentage of maximum used to get ToT level
            TOT30_begin = 0, TOT30_end = 0;

            getTOTs( fast_waveform, TOT20_div, TOT30_div,
                TOT20_begin, TOT20_end, TOT30_begin, TOT30_end,
                max_value_a, time_fast_wf);

            TOT20.push_back( (TOT20_end-TOT20_begin) );           
            TOT30.push_back( (TOT30_end-TOT30_begin) );

            
            //-----------  Getting information for BAT  -------------------------------------------------------------------------//


            sliceWaveform_BAT(fast_waveform, integrals_slices, matching_slices, TOT20_begin, TOT20_end);
            

            //-----------  Getting travelled Z  -----------//
            //    **SHOULD CORRECT FOR MINIMUM TILTNESS OFFSET
            

            delta_z = (TOT20_end - TOT20_begin) * (4/3) * (drift_vel/1000.); // in cm,         
            travelled_Z.push_back(delta_z);


            //-----------  Determination of track quadrant from the PMT   -------------------------------------------------------//


            wf_integral = accumulate(fast_waveform->begin(), fast_waveform->end(), 0);
            integrals_wfs.push_back(wf_integral);


            //----------- Calculation of skewness of the Bragg Peak  ------------------------------------------------------------//

            wf_peaks.clear();
            getSkewness_BraggPeak(fast_waveform, time_fast_wf, TOT20_begin, TOT20_end, 
            max_value_a, max_value_t, skew_ratio, wf_peaks);


            //----------- Make final waveform plot with all the information  ---------------------------------------------------//


            create_and_print_wf_graph_lines2(Form("WF_run_%i_evt_%i_trg_%i_ch_%i",pmt_run,pmt_event,pmt_trigger,pmt_channel), 
            time_fast_wf, fast_waveform,
            TOT20_begin, TOT20_end, max_value_a * TOT20_div, 
            TOT30_begin, TOT30_end, max_value_a * TOT30_div, 
            wf_peaks);


            //----------- Calculation of all the final variables of interest  --------------------------------------------------//
           
            if ( pmt_channel == 4 ) {                   // Perform actions only when I have the 4 waveforms
                
                cout << "\nPMT Track information: \n" << endl; 


                //----------- Calculate direction in Z, and respective score  ------------------------//

                direction = 0, dir_score = 0;
                getDirectionScore(skew_ratio, direction, dir_score, false);

                if      (direction == -1 ) cout << "--> Moving towards the GEMs with score: "     << dir_score*100.   << endl;
                else if (direction == +1 ) cout << "--> Moving towards the cathode with score: "  << dir_score*100.   << endl;
                else if (direction ==  0 ) cout << "--> Ambiguous. Score: " << dir_score*100. << endl;
            

                //----------- Calculate travelled Z  ------------------------//

                avg_travel_z = accumulate(travelled_Z.begin(), travelled_Z.end(), 0.0) / travelled_Z.size();
                cout << "--> The average travelled Z (cm) is: " << avg_travel_z << endl;


                //----------- Calculate Quadrant ---------------------------//
                
                getQuadrantPMT(integrals_wfs, quadrant_pmt);
                cout << "--> The track is in the quadrant: " << quadrant_pmt << endl;


                //----------- Particle ID ---------------------------------//

                pmt_PID_total = false;
                getAlphaIdentification(TOT20, TOT30, pmt_PID_total, true);
                cout << "--> The particle in this trigger was identified as an alpha: " << pmt_PID_total << endl;


                //----------- Run BAT routine ----------------------------//

                create_bat_input( pmt_run, pmt_event, pmt_trigger, integrals_slices, "input_for_bat.txt");
                run_bat("input_for_bat.txt", "output_from_bat.txt");
                read_bat("output_from_bat.txt",points_bat, true);
                print_BAT_CAM_match(h_fullsize, points_bat, "bat", "BAT Match");



                //----------- Collect all the relevant info for posterior analysis  -----------//

                PMT_alphas.push_back({
                    .run = pmt_run,
                    .pic  = pmt_event,
                    .trg = pmt_trigger,

                    .dir = direction,
                    .prob = (dir_score*100.),
                    .trv_Z = avg_travel_z,

                    .quad = quadrant_pmt,

                    .its_alpha = pmt_PID_total
                });


                //----------- Clear vector used for the 4 waveforms  -----------//

                TOT30.clear();
                TOT20.clear();                
                skew_ratio.clear();                
                integrals_wfs.clear();                
                wf_peaks.clear();                
                travelled_Z.clear();
                integrals_slices.clear();
                points_bat.clear();
                points_cam.clear();
            }
        }
    }
    reco_data_pmt->Close();

    /* **************************************  COMBINED ANALYSIS  ********************************************************************  */

    /* ************* Combined information ********** */

    int run, picture, trigger;
    double IP_X, IP_Y, IP_Z;
    double track_end_X, track_end_Y, track_end_Z;
    double Z_angle, XY_angle;
    double Z_length, XY_length;
    double full_length;
    double pmt_direction;
    double pmt_direction_score;
    int cam_quad, pmt_quad;
    // bool cam_PID, pmt_PID;
    double cam_energy, pmt_energy;

    TTree *tree_3D = new TTree("3D_Alpha_Tracks", "3D Alpha Tracks");

    tree_3D->Branch("run", &run, "run/I");
    tree_3D->Branch("picture", &picture, "picture/I");
    tree_3D->Branch("trigger", &trigger, "trigger/I");
    tree_3D->Branch("IP_X", &IP_X, "IP_X/D");
    tree_3D->Branch("IP_Y", &IP_Y, "IP_Y/D");
    tree_3D->Branch("IP_Z", &IP_Z, "IP_Z/D");
    tree_3D->Branch("track_end_X", &track_end_X, "track_end_X/D");
    tree_3D->Branch("track_end_Y", &track_end_Y, "track_end_Y/D");
    tree_3D->Branch("track_end_Z", &track_end_Z, "track_end_Z/D");
    tree_3D->Branch("Z_angle", &Z_angle, "Z_angle/D");
    tree_3D->Branch("XY_angle", &XY_angle, "XY_angle/D");
    tree_3D->Branch("Z_length", &Z_length, "Z_length/D");
    tree_3D->Branch("XY_length", &XY_length, "XY_length/D");
    tree_3D->Branch("full_length", &full_length, "full_length/D");
    tree_3D->Branch("pmt_direction", &pmt_direction, "pmt_direction/D");
    tree_3D->Branch("pmt_direction_score", &pmt_direction_score, "pmt_direction_score/D");
    tree_3D->Branch("cam_quad", &cam_quad, "cam_quad/I");
    tree_3D->Branch("pmt_quad", &pmt_quad, "pmt_quad/I");
    // tree_3D->Branch("cam_PID", &cam_PID, "cam_PID/O");
    // tree_3D->Branch("pmt_PID", &pmt_PID, "pmt_PID/O");
    tree_3D->Branch("cam_energy", &cam_energy, "cam_energy/D");
    tree_3D->Branch("pmt_energy", &pmt_energy, "pmt_energy/D");

    for (const auto& cam : CAM_alphas) {

        for (const auto& pmt : PMT_alphas) {
        
            if ( pmt.run == cam.run ) {

                if ( pmt.pic == cam.pic ) {
                
                    if (pmt.quad == cam.quad ) {        

                        if ( pmt.its_alpha == true && cam.its_alpha == true) {

                            cout << "\n# Matched alpha in quadrant: " << pmt.quad << "; in trigger: " << pmt.trg << "; with Alpha-PID = " << pmt.its_alpha << endl;
                            
                            //----------- Creating 3D alpha track information  -----------//

                            IP_X = cam.IP_X_cm;
                            IP_Y = cam.IP_Y_cm;               
                            IP_Z = 25.0;                          // Absolute Z fixed.

                            track_end_X = IP_X + ( cam.trv_XY * cos(cam.angle_XY * TMath::Pi()/180.));
                            track_end_Y = IP_Y + ( cam.trv_XY * sin(cam.angle_XY * TMath::Pi()/180.));
                            track_end_Z = IP_Z + ( pmt.trv_Z  * pmt.dir);

                            Z_angle = atan(pmt.trv_Z/cam.trv_XY) * 180. / TMath::Pi();

                            full_length = TMath::Sqrt(pow(cam.trv_XY,2) + pow(pmt.trv_Z,2));

                            cout << "\n\t ** 3D Alpha track information: ** \n" << endl; 
                            cout << "--> Position, X: " << IP_X << "; Y: " << IP_Y  << endl;
                            cout << "--> Travelled XY: "    << cam.trv_XY   << endl;
                            cout << "--> Angle XY (#phi): "        << cam.angle_XY << endl;
                            cout << "--> Travelled Z: "     << pmt.trv_Z    << endl;
                            cout << "--> Direction in Z: "  << pmt.dir      << " at " << abs(pmt.prob) << " score" << endl;
                            cout << "--> Angle Z (#theta): "     << Z_angle  << endl;
                            cout << "--> 3D alpha length (cm): " << full_length  << endl;

                            file_root->cd(Form("Run_%i_ev_%i", cam.run, cam.pic));
                            build_3D_vector(IP_X,track_end_X,IP_Y,track_end_Y,IP_Z,track_end_Z,cam.trv_XY, cam.angle_XY, pmt.trv_Z, pmt.dir, pmt.prob, Z_angle, full_length, pmt.run, pmt.pic, pmt.trg);
                        
                            tree_3D->Fill();

                        } else continue;
                    } else continue;
                } else continue;
            } else continue;
        } //close for pmt
    } //close for cam

    file_root->cd();
    tree_3D->Write();

    file_root->Close();
    cout << "**Finished**" << endl;
    myapp->Run();

    return 0;
}


/*  ***************************************************************************************************************************  */
/*  ***************************************************************************************************************************  */
/*  ***************************************************************************************************************************  */


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