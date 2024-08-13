#include <iostream>
#include <numeric>
#include "TFile.h"
// #include "TH1F.h"
#include "TH2F.h"
// #include "TF1.h"
// #include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TApplication.h"
#include "TLine.h"
#include "TTree.h"
#include "TROOT.h"
#include "TFitResult.h"

#include "TKey.h"

#include "../Track_analyzer/Analyzer.h"   //** Change the path to the correct one
#include "./include.h"

using namespace std;

void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E);
void deleteNonAlphaDirectories(const char *filename, bool deleteAll = false);
std::string exec(const char* cmd); 

int main(int argc, char**argv) {
    auto start = std::chrono::high_resolution_clock::now();

    /* ****************************************  Running options   **************************************************************************  */

    bool save_everything = false;    //opposite of saving *only* the alpha-tree

    string mode = argv[ 1 ]; // debug or full
    bool batch_mode;
    if (mode == "full") batch_mode = true;
    if (mode == "debug") batch_mode = false, save_everything = true;

    /* ****************************************  Inputs and outputs   **************************************************************************  */

    if(batch_mode) gROOT->SetBatch(kTRUE);

    if (argc < 4) {
        cout << "Usage: ./Alpha_3D_LIME mode <CAMERA reco file> <PMT reco file> <output file> " << endl;
        return 1;
    }

    std::cout << std::boolalpha; // Enable boolalpha
    TApplication *myapp=new TApplication("myapp",0,0); cout << "\n" << endl;
    
    string filename_cam = argv[ 2 ];

    string filename_pmt = argv[ 3 ];

    string outputfile = argv[ 4 ];
    string final_out; 
    if (mode == "full") final_out = "out/" + outputfile + ".root";
    if (mode == "debug") final_out = "out/debug_" + outputfile + ".root";
    TFile* file_root = new TFile(final_out.c_str(),"recreate");

    //for debugging and testing
    int debug_event = 0;
    if (mode == "debug") debug_event = atof(argv[ 5 ]);
    // int debug_trigger = atof(argv[ 5 ]);

    /* ****************************************  Definition of variables   **************************************************************************  */

    // -----------------  Detector / Physics constants  ----------------- //

    const double drift_vel = 5.471; // cm/µs // electron drift velocity in LIME, He:CF4 60:40, 800 V/cm 
    const double granularity = 0.0155; // cm/pixel 

    // -----------------  Analyzer constants  ----------------- //

    // For alphas (David)
    Int_t NPIP=2000; Float_t wFac=3.5;
    
    // For ER (original from Samuele)
    // Int_t NPIP=80; Float_t wFac=2.5;

    
    // Cluster pixels reading
    int nSc_red=0;
    vector<int> BeginScPix;
    vector<int> EndScPix;

    // X axis for the waveforms
    vector<int> time_fast_wf(1024); iota(time_fast_wf.begin(), time_fast_wf.end(), 0);
    vector<int> time_slow_wf(4000); iota(time_slow_wf.begin(), time_slow_wf.end(), 0);

    // Direction
    int direction; // -1 = towards GEM ; 1 = towards cathode; 0 = ambiguous
    double dir_score;
    bool verbose_dir_score = false;
    
    //Waveform analysis
    double wf_integral;
    vector<double> integrals_wfs;
    double avg_skew;
    vector<double> skew_ratio;          // ** check how to calculate skewness. Useful to then connect with Bayes fit.
    vector<pair<int,double>> wf_peaks;
    vector<pair<int,double>> wf_peaks_energy_dep;
    double wf_Npeaks_ed = 0;
    int matching_slices = 5;
    vector<vector<double>> integrals_slices;

    // Travelled Z
    double delta_z;
    double avg_travel_z;
    vector<double> travelled_Z;

    // Variables for ToT
    int TOT20_begin = 0, TOT20_end = 0; double TOT20_div; vector<double> TOT20;
    int TOT30_begin = 0, TOT30_end = 0; double TOT30_div; vector<double> TOT30;

    // PID
    int quadrant_pmt;
    int quadrant_cam;
    bool pmt_PID_total;
    bool cam_PID = false;

    // BAT and CAM matching
    vector<pair<double,double>> points_bat;
    vector<pair<double,double>> points_cam;

    // CAM variables
    double x_impact, y_impact;
    double xbar, ybar;
    double angle_cam;

    // Profiles
    double fitAmp, fitMean, fitSigma, fitConst;
    double fitAmpError, fitMeanError, fitSigmaError, fitConstError;

    // PMT
    double fitted_lum;

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
    vector<float> sc_xmin;        sc_xmin.reserve(150000);        tree_cam->SetBranchAddress("sc_xmin",       sc_xmin.data());  
    vector<float> sc_ymin;        sc_ymin.reserve(150000);        tree_cam->SetBranchAddress("sc_ymin",       sc_ymin.data());  
    vector<float> sc_nhits;       sc_nhits.reserve(150000);       tree_cam->SetBranchAddress("sc_nhits",      sc_nhits.data());  

    UInt_t nSc;                     tree_cam->SetBranchAddress("nSc",   &nSc);
    UInt_t Nredpix=0;               tree_cam->SetBranchAddress("nRedpix",&Nredpix);

    // --> Reserve even higher than the expected number of pixels to avoid problems with the vector size
    vector<float> sc_redpixID;      sc_redpixID.reserve(1500000);    tree_cam->SetBranchAddress("sc_redpixIdx",sc_redpixID.data());
    // Run 3
    // vector<int> XPix;               XPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_iy",       XPix.data());  
    // vector<int> YPix;               YPix.reserve(150000);           tree_cam->SetBranchAddress("redpix_ix",       YPix.data());
    // Run 4
    vector<int> XPix;               XPix.reserve(1500000);           tree_cam->SetBranchAddress("redpix_ix",       XPix.data());  
    vector<int> YPix;               YPix.reserve(1500000);           tree_cam->SetBranchAddress("redpix_iy",       YPix.data());
    
    vector<float> ZPix;             ZPix.reserve(1500000);           tree_cam->SetBranchAddress("redpix_iz",       ZPix.data());  


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
    
    cout << "\n\n************   Analysis  CAMERA    ************\n\n" << endl;

    vector<AlphaTrackCAM> CAM_alphas;
    Int_t nentries = (Int_t)tree_cam->GetEntries();
    for (Int_t cam_entry = 0; cam_entry < nentries; ++cam_entry) {
        
        tree_cam->GetEntry(cam_entry);

        if (mode == "debug" && cam_event != debug_event) continue;

        sc_redpixID.clear();
        ScIndicesElem(nSc, Nredpix, sc_redpixID.data(), nSc_red, BeginScPix, EndScPix);

        for ( int sc_i = 0; sc_i < nSc_red; sc_i++ ) {

            cam_PID = false;
            cout << "\n\n\t==> Cam run: " << cam_run << "; event: " << cam_event << "; cluster ID: " << sc_i << "\n\n" << endl;

            // if ( (sc_rms[nc] > 6) && (0.152 * sc_tgausssigma[nc] > 0.3) && (0.152 * sc_length[nc] < 80) && (sc_integral[nc] > 1000) && (sc_width[nc] / sc_length[nc] > 0.8)            //cam cut
            // && ( ( (sc_xmean[nc]-1152)*(sc_xmean[nc]-1152) + (sc_ymean[nc]-1152)*(sc_ymean[nc]-1152) ) <(800*800)) ) {
            if ( sc_integral[sc_i]/sc_nhits[sc_i] > 25 && sc_length[sc_i] > 100 && sc_width[sc_i] > 50 ) {   //Alpha cut fropm Giorgio
            // if ( sc_integral[sc_i]/sc_nhits[sc_i] < 25 && sc_length[sc_i] > 100 && sc_width[sc_i] > 0 ) {   //Cosmics cut! (for thesis)
            // if ( sc_integral[sc_i]/sc_nhits[sc_i] < 25 && sc_length[sc_i] > 10 && sc_width[sc_i] > 0 ) {   //Cosmics cut! (for thesis)

                if (save_everything) {
                    TString folderName = Form("Run_%i_ev_%i", cam_run, cam_event);
                    if (!file_root->GetDirectory(folderName)) file_root->mkdir(folderName);
                    file_root->cd(Form("Run_%i_ev_%i", cam_run, cam_event));
                }

                cam_PID = true;

                //----------- Build Analyser track  -----------//

                Analyzer Track(Form("Track_run_%i_ev_%i_cl_%i", cam_run, cam_event,sc_i),XPix.data(),YPix.data(),ZPix.data(),BeginScPix[sc_i],EndScPix[sc_i]);
                
                Track.SetWScal(wFac);
                Track.SetNPIP(NPIP);
                Track.ApplyThr();
                // // Track.RemoveNoise(25); //iron
                // Track.RemoveNoise(25); //cosmics
                // Track.RemoveNoise(75);
                Track.RemoveNoise(100); //alphas
                Track.ImpactPoint(Form("IP_run_%i_ev_%i_cl_%i", cam_run, cam_event,sc_i));
                Track.ScaledTrack(Form("Rebinned_run_%i_ev_%i_cl_%i", cam_run, cam_event,sc_i));
                Track.Direction();
                Track.ImprCorrectAngle();
                Track.BuildLineDirection();

                points_cam = Track.GetLinePoints(matching_slices,"edges");

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

                if(save_everything) {
                    const char* name = Form("ev_%i_run_%i_cluster_%i", cam_run, cam_event,sc_i);
                    Track.PlotandSavetoFileCOLZ_fullSize(Form("Track_%s",name));     
                    Track.PlotandSavetoFileDirectionalFull(Form("X_Y_Analyser_%s",name));
                    
                    TFitResultPtr fitResult; 
                    printTrackProfilesAndFit(Track.FillProfile(false,name), Track.FillProfile(true,name),Form("Profiles_%s",name), fitResult);     //"memory leak" because two TH2D are created with the same name.
                    fitAmp   = fitResult->Parameter(0), fitAmpError   = fitResult->ParError(0); 
                    fitMean  = fitResult->Parameter(1), fitMeanError  = fitResult->ParError(1); 
                    fitSigma = fitResult->Parameter(2), fitSigmaError = fitResult->ParError(2);
                    fitConst = fitResult->Parameter(3), fitConstError = fitResult->ParError(3);
                }

                //----------- Collect all the relevant info for posterior analysis  -----------//

                CAM_alphas.push_back({
                    .run = cam_run,
                    .pic  = cam_event,
                    .cluster = sc_i,

                    .angle_XY = angle_cam,
                    .trv_XY = sc_length[sc_i] * granularity,

                    .quad = quadrant_cam,

                    .IP_X_cm = x_impact,
                    .IP_Y_cm = y_impact,

                    .its_alpha = cam_PID,

                    .track_cam = points_cam,

                    .energy      = sc_integral[sc_i],
                    .nhits       = sc_nhits[sc_i],
                    .width       = sc_width[sc_i] * granularity,
                    .xmean       = sc_xmean[sc_i] * granularity,
                    .ymean       = sc_ymean[sc_i] * granularity,
                    .rms         = sc_rms[sc_i],
                    .tgausssigma = sc_tgausssigma[sc_i],

                    .fitSig = fitSigma

                });      

                //----------- Cleanup ----------------- -----------//

                points_cam.clear();

            } else cout << "--> The particle in this cluster was identified as an alpha: " << cam_PID << endl;
        }
        sc_redpixID.resize(nSc); // Resize to avoid problem with vector sizes
    }
    reco_data_cam->Close();
    

    /* *********************************************  Analysis PMT  *******************************************  */
    
    cout << "\n\n************   Analysis  PMT    ************\n\n" << endl;

    vector<AlphaTrackPMT> PMT_alphas;
    Int_t nentries_pmt = (Int_t)tree_pmt->GetEntries();
    for (Int_t pmt_wf = 0; pmt_wf < nentries_pmt; pmt_wf++) {  // Each entry is one waveform
        
        tree_pmt->GetEntry(pmt_wf);

        if ( pmt_sampling != 1024)                                    continue; 
        if ( mode == "debug" && pmt_event != debug_event )            continue;
        if ( !found_clusters_in_evt(CAM_alphas, pmt_run, pmt_event) ) continue;

        if ( pmt_channel == 1) cout << "\n\n\t==> PMT run: " << pmt_run << "; event: " << pmt_event << "; trigger: " << pmt_trigger << "; sampling: " << pmt_sampling << endl;
        if (save_everything) file_root->cd(Form("Run_%i_ev_%i", pmt_run, pmt_event));


        //-----------  Put branch information into vector for further analysis  --------------------------------------------//

        shared_ptr<vector<double>> fast_waveform = make_shared<vector<double>>();
        for (int j = 0; j < maxArraySize_fast; ++j) fast_waveform->push_back(pmt_waveforms[j]);
        // shared_ptr<vector<double>> slow_waveform = make_shared<vector<double>>();


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

        //----------- Calculation of nPeaks for PID  ------------------------------------------------------------------------//

        wf_peaks_energy_dep.clear();
        findPeaks(fast_waveform, 50 , wf_peaks_energy_dep);
        wf_Npeaks_ed += wf_peaks_energy_dep.size();


        //----------- Make final waveform plot with all the information  ---------------------------------------------------//

        if(save_everything){ 

            create_and_print_wf_graph_lines3(Form("WF_run_%i_evt_%i_trg_%i_ch_%i",pmt_run,pmt_event,pmt_trigger,pmt_channel), 
            time_fast_wf, fast_waveform,
            TOT20_begin, TOT20_end, max_value_a * TOT20_div, 
            TOT30_begin, TOT30_end, max_value_a * TOT30_div, 
            wf_peaks, wf_peaks_energy_dep);
        }

        //----------- Calculation of all the final variables of interest  --------------------------------------------------//
        
        if ( pmt_channel == 4 ) {                   

            cout << "\n\n**PMT Track information: \n" << endl; 

            //----------- Particle ID ---------------------------------//

            pmt_PID_total = false;
            getAlphaIdentification(TOT20, TOT30, wf_Npeaks_ed, pmt_PID_total, true);
            cout << "--> The particle in this trigger was identified as an alpha: " << pmt_PID_total << endl;

            
            if (pmt_PID_total) {

                //----------- Calculate direction in Z, and respective score  ------------------------//

                direction = 0, dir_score = 0;


                if (!skew_ratio.empty()) {
                    getDirectionScore(skew_ratio, direction, dir_score, false);
                } else {
                    direction = 0; dir_score = 0;
                    cout << "No skewness calculated. Direction is ambiguous." << endl;
                }

                if      (direction == -1 ) cout << "--> Moving towards the GEMs with score: "     << dir_score*100.   << endl;
                else if (direction == +1 ) cout << "--> Moving towards the cathode with score: "  << dir_score*100.   << endl;
                else if (direction ==  0 ) cout << "--> Ambiguous. Score: " << dir_score*100. << endl;


                //----------- Calculate travelled Z  ------------------------//

                avg_travel_z = accumulate(travelled_Z.begin(), travelled_Z.end(), 0.0) / travelled_Z.size();
                cout << "--> The average travelled Z (cm) is: " << avg_travel_z << endl;


                //----------- Calculate Quadrant ---------------------------//
                
                getQuadrantPMT(integrals_wfs, quadrant_pmt);
                cout << "--> The track is in the quadrant: " << quadrant_pmt << endl;


                //----------- Run BAT routine ----------------------------//

                fitted_lum = 0;
                /* File naming a bit hardcoded here */
                create_bat_input( pmt_run, pmt_event, pmt_trigger, integrals_slices, "bat_files/input_for_bat.txt");
                run_bat("bat_files/input_for_bat.txt", "bat_files/output_from_bat.txt", "../BAT_PMTs/./runfit.out");
                read_bat("bat_files/output_from_bat.txt",points_bat, fitted_lum, false);

                //----------- Collect all the relevant info for posterior analysis  -----------//

                PMT_alphas.push_back({
                    .run = pmt_run,
                    .pic  = pmt_event,
                    .trg = pmt_trigger,

                    .dir = direction,
                    .prob = (dir_score*100.),
                    .trv_Z = avg_travel_z,

                    .quad = quadrant_pmt,

                    .its_alpha = pmt_PID_total,

                    .track_pmt = points_bat,

                    .energy = fitted_lum,

                    .num_peaks = wf_Npeaks_ed/4 
                });
            }

            //----------- Clear vector used for the 4 waveforms  -----------//

            TOT30.clear();
            TOT20.clear();                
            skew_ratio.clear();                
            integrals_wfs.clear();                
            wf_peaks.clear();                
            wf_peaks_energy_dep.clear();                
            travelled_Z.clear();
            integrals_slices.clear();
            points_bat.clear();
            wf_Npeaks_ed = 0;
        }
    }
    reco_data_pmt->Close();

    /* **************************************  COMBINED ANALYSIS  ********************************************************************  */

    cout << "\n\n************   COMBINED Analysis   ************\n\n" << endl;

    /* ************* Combined information ********** */

    // General
    int run, picture, trigger;
    
    //3D track variables
    double IP_X, IP_Y, IP_Z;
    double track_end_X, track_end_Y, track_end_Z;
    double Z_angle, XY_angle;
    double Z_length, XY_length;
    double full_length;
    
    // Head-tail
    double pmt_direction;
    double pmt_direction_score;
    
    //Association
    int cam_quad, pmt_quad;
    bool cam_ParticleID, pmt_ParticleID;
    vector<pair<double,double>> distances;
    vector<tuple<double, size_t, size_t>> small_dist_assoc;

    // Remaining variables
    double pmt_energy;
    double pmt_peaks;
    double cam_energy;
    double cam_nhits;
    double cam_width;
    double cam_xmean;
    double cam_ymean;
    double cam_rms;
    double cam_tgausssigma;

    TTree *tree_3D = new TTree("AlphaEvents", "3D Alpha Tracks");

    // General
    tree_3D->Branch("run", &run, "run/I");
    tree_3D->Branch("picture", &picture, "picture/I");
    tree_3D->Branch("trigger", &trigger, "trigger/I");

    // 3D track variables
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

    // Head-tail
    tree_3D->Branch("pmt_direction", &pmt_direction, "pmt_direction/D");
    tree_3D->Branch("pmt_direction_score", &pmt_direction_score, "pmt_direction_score/D");

    // Association
    tree_3D->Branch("cam_quad", &cam_quad, "cam_quad/I");
    tree_3D->Branch("pmt_quad", &pmt_quad, "pmt_quad/I");
    tree_3D->Branch("cam_ParticleID", &cam_ParticleID, "cam_ParticleID/O");
    tree_3D->Branch("pmt_ParticleID", &pmt_ParticleID, "pmt_ParticleID/O");
    tree_3D->Branch("distances", &distances);

    // Remaining variables
    tree_3D->Branch("pmt_energy", &pmt_energy, "pmt_energy/D");
    tree_3D->Branch("pmt_peaks", &pmt_peaks, "pmt_peaks/D");
    tree_3D->Branch("cam_energy", &cam_energy, "cam_energy/D");
    tree_3D->Branch("cam_nhits", &cam_nhits, "cam_nhits/D");
    tree_3D->Branch("cam_width", &cam_width, "cam_width/D");
    tree_3D->Branch("cam_xmean", &cam_xmean, "cam_xmean/D");
    tree_3D->Branch("cam_ymean", &cam_ymean, "cam_ymean/D");
    tree_3D->Branch("cam_rms", &cam_rms, "cam_rms/D");
    tree_3D->Branch("cam_tgausssigma", &cam_tgausssigma, "cam_tgausssigma/D");
    
    // Capture the current Git commit hash
    string git_commit_hash = exec("git rev-parse HEAD");

    // Create a branch in the tree to store the commit hash
    char git_commit_hash_cstr[41]; // Git commit hash is 40 characters long
    strncpy(git_commit_hash_cstr, git_commit_hash.c_str(), 40);
    git_commit_hash_cstr[40] = '\0'; // Ensure null-termination
    tree_3D->Branch("git_commit_hash", git_commit_hash_cstr, "git_commit_hash/C");

///////////////////////////////////

    // To make the BAT association, I need to re-arrange some data because I'm saving the alphas individually as "entries"
    // So I now re-group the events per picture number in a map so that I can ASSOCIATE trigger and clusters in the same event.  

    // Map where the key will be the {run,pic} number and the value the alpha(s) in that pic.
    std::map<RunPicKey, std::vector<decltype(CAM_alphas)::value_type>> cam_map;
    std::map<RunPicKey, std::vector<decltype(PMT_alphas)::value_type>> pmt_map;

    // Group alphas by run number and picture number
    for (const auto& cam : CAM_alphas) {
        RunPicKey key = {cam.run, cam.pic};
        cam_map[key].push_back(cam);
    }
    for (const auto& pmt : PMT_alphas) {
        RunPicKey key = {pmt.run, pmt.pic};
        pmt_map[key].push_back(pmt);
    }

    cout << "=> Total amount of alphas in the CAMERA: " << CAM_alphas.size() << endl;
    cout << "=> Total amount of alphas in the PMT: "    << PMT_alphas.size() << endl;
    cout << endl;

    for (const auto& cam_entry : cam_map) {

        const RunPicKey& key = cam_entry.first;
        const auto& cam_entries = cam_entry.second;
        const auto& pmt_entries = pmt_map[key]; // Corresponding PMT entries for the current run number and picture number

        small_dist_assoc.clear();

        one_to_one_association(small_dist_assoc, pmt_entries, cam_entries, false);

        cout << "Associations for run " << key.run << ", picture " << key.pic << endl;
        cout << "# CAM alphas: " << cam_entries.size() << ";  # PMT alphas: " << pmt_entries.size() << endl;

        // cout << "Total # combinatorics: " << small_dist_assoc.size() << endl;
        // for (size_t i = 0; i < small_dist_assoc.size(); ++i) {
        //     double dd = std::get<0>(small_dist_assoc[i]);
        //     size_t cam_idx = std::get<1>(small_dist_assoc[i]);
        //     size_t pmt_idx = std::get<2>(small_dist_assoc[i]);
        //     const auto& cam = cam_entries[cam_idx];
        //     const auto& pmt = pmt_entries[pmt_idx];
        //     cout << "Possible combinations: " << endl;
        //     cout << " -> CAM cluster # = " << cam_idx << "; PMT trigger" << pmt.trg << "; Distance: " << dd << endl;
        // }

        int n_assoc = min(cam_entries.size(), pmt_entries.size());
        cout << "==> We'll perform > " << n_assoc << " < associations." << endl;
        
        for (size_t assoc_i = 0; assoc_i < n_assoc; ++assoc_i) {
            
            double dd = std::get<0>(small_dist_assoc[assoc_i]);
            size_t cam_idx = std::get<1>(small_dist_assoc[assoc_i]);
            size_t pmt_idx = std::get<2>(small_dist_assoc[assoc_i]);
            
            const auto& cam = cam_entries[cam_idx];
            const auto& pmt = pmt_entries[pmt_idx];
        
            cout << endl;
            cout << "---------------------------------------------------\n" << endl;
            cout << "*** Association # "    << assoc_i+1<< ":"              << endl;
            cout << endl;
            cout << " -> CAM cluster = "    << cam.cluster  << endl;
            cout << " -> PMT trigger "      << pmt.trg      << endl;
            cout << " -> Distance: "        << dd           << endl;

            run = cam.run; picture = cam.pic; trigger = pmt.trg;
            cout << Form("\n ==> Event in run %i, event %i, trigger %i, quadrant %i, with Alpha-PID = %i", cam.run, cam.pic, pmt.trg, pmt.quad, pmt.its_alpha) << endl;
            
            //----------- Creating 3D alpha track information  -----------//

            IP_X = cam.IP_X_cm;
            IP_Y = cam.IP_Y_cm;               
            IP_Z = 25.0;                          // Absolute Z fixed.

            track_end_X = IP_X + ( cam.trv_XY * cos(cam.angle_XY * TMath::Pi()/180.));
            track_end_Y = IP_Y + ( cam.trv_XY * sin(cam.angle_XY * TMath::Pi()/180.));
            track_end_Z = IP_Z + ( pmt.trv_Z  * pmt.dir);                               // if dir = 0, track doesn't show Z direction

            Z_angle = atan(pmt.trv_Z/cam.trv_XY) * 180. / TMath::Pi();
            XY_angle = cam.angle_XY;

            Z_length = pmt.trv_Z;
            XY_length = cam.trv_XY;
            full_length = TMath::Sqrt(pow(cam.trv_XY,2) + pow(pmt.trv_Z,2));

            pmt_direction = pmt.dir;
            pmt_direction_score = pmt.prob;

            cam_quad = cam.quad;
            pmt_quad = pmt.quad;

            cam_ParticleID = cam.its_alpha;
            pmt_ParticleID = pmt.its_alpha;


            //-----------  Accuracy of BAT  -------------------------------//

            calculate_distance(pmt.track_pmt, cam.track_cam, distances, false);
            
            //----------- Verbose information  -----------//

            cout << "\n  ** 3D Alpha track information: ** \n" << endl; 
            cout << "--> Position, X: " << IP_X << "; Y: " << IP_Y  << endl;
            cout << "--> Travelled XY: "    << XY_length   << endl;
            cout << "--> Angle XY (#phi): "        << XY_angle << endl;
            cout << "--> Travelled Z: "     << Z_length    << endl;
            cout << "--> Direction in Z: "  << pmt_direction      << " at " << abs(pmt_direction_score) << " score" << endl;
            cout << "--> Angle Z (#theta): "     << Z_angle  << endl;
            cout << "--> 3D alpha length (cm): " << full_length  << endl;
            cout << endl;
            cout << "---------------------------------------------------" << endl;

            if(save_everything){
                file_root->cd(Form("Run_%i_ev_%i", cam.run, cam.pic));
                Points_BAT_CAM(pmt.track_pmt, cam.track_cam, Form("BAT_CAM_association_#%i", assoc_i));
                build_3D_vector(IP_X,track_end_X,IP_Y,track_end_Y,IP_Z,track_end_Z,cam.trv_XY, cam.angle_XY, pmt.trv_Z, pmt.dir, abs(pmt_direction_score), Z_angle, full_length, pmt.run, pmt.pic, pmt.trg, true, cam.fitSig*granularity);
            }
            //-----------  Variables for mixed analysis  -------------------------------//

            pmt_energy = pmt.energy;
            pmt_peaks = pmt.num_peaks;
            cam_energy = cam.energy;
            cam_nhits = cam.nhits;
            cam_width = cam.width;
            cam_xmean = cam.xmean;
            cam_ymean = cam.ymean;
            cam_rms = cam.rms;
            cam_tgausssigma = cam.tgausssigma;

            //-----------  Tree filling and variables cleaning  ----------//
        
            tree_3D->Fill();
            distances.clear();

        }

    }

    /*
    // for (const auto& cam : CAM_alphas) {

        // for (const auto& pmt : PMT_alphas) {
        
            if ( pmt.run == cam.run && pmt.pic == cam.pic) {

                if (pmt.quad == cam.quad ) {         //quad selection is too sensitive to alphas in the center    
                    if ( pmt.its_alpha == true && cam.its_alpha == true) {

                        run = cam.run; picture = cam.pic; trigger = pmt.trg;

                        cout << Form("\n\n ==> Matched alpha in run %i, event %i, trigger %i, quadrant %i, with Alpha-PID = %i", cam.run, cam.pic, pmt.trg, pmt.quad, pmt.its_alpha) << endl;
                        
                        //----------- Creating 3D alpha track information  -----------//

                        IP_X = cam.IP_X_cm;
                        IP_Y = cam.IP_Y_cm;               
                        IP_Z = 25.0;                          // Absolute Z fixed.

                        track_end_X = IP_X + ( cam.trv_XY * cos(cam.angle_XY * TMath::Pi()/180.));
                        track_end_Y = IP_Y + ( cam.trv_XY * sin(cam.angle_XY * TMath::Pi()/180.));
                        track_end_Z = IP_Z + ( pmt.trv_Z  * pmt.dir);                               // if dir = 0, track doesn't show Z direction

                        Z_angle = atan(pmt.trv_Z/cam.trv_XY) * 180. / TMath::Pi();
                        XY_angle = cam.angle_XY;

                        Z_length = pmt.trv_Z;
                        XY_length = cam.trv_XY;
                        full_length = TMath::Sqrt(pow(cam.trv_XY,2) + pow(pmt.trv_Z,2));

                        pmt_direction = pmt.dir;
                        pmt_direction_score = pmt.prob;

                        cam_quad = cam.quad;
                        pmt_quad = pmt.quad;

                        cam_ParticleID = cam.its_alpha;
                        pmt_ParticleID = pmt.its_alpha;


                        //-----------  Accuracy of BAT  -------------------------------//

                        calculate_distance(pmt.track_pmt, cam.track_cam, distances, true);
                        
                        //----------- Verbose information  -----------//

                        cout << "\n\t ** 3D Alpha track information: ** \n" << endl; 
                        cout << "--> Position, X: " << IP_X << "; Y: " << IP_Y  << endl;
                        cout << "--> Travelled XY: "    << XY_length   << endl;
                        cout << "--> Angle XY (#phi): "        << XY_angle << endl;
                        cout << "--> Travelled Z: "     << Z_length    << endl;
                        cout << "--> Direction in Z: "  << pmt_direction      << " at " << abs(pmt_direction_score) << " score" << endl;
                        cout << "--> Angle Z (#theta): "     << Z_angle  << endl;
                        cout << "--> 3D alpha length (cm): " << full_length  << endl;

                        if(save_everything){
                            file_root->cd(Form("Run_%i_ev_%i", cam.run, cam.pic));
                            build_3D_vector(IP_X,track_end_X,IP_Y,track_end_Y,IP_Z,track_end_Z,cam.trv_XY, cam.angle_XY, pmt.trv_Z, pmt.dir, pmt.prob, Z_angle, full_length, pmt.run, pmt.pic, pmt.trg);
                        }
                        //-----------  Variables for mixed analysis  -------------------------------//

                        pmt_energy = pmt.energy;
                        cam_energy = cam.energy;
                        cam_nhits = cam.nhits;
                        cam_width = cam.width;
                        cam_xmean = cam.xmean;
                        cam_ymean = cam.ymean;
                        cam_rms = cam.rms;
                        cam_tgausssigma = cam.tgausssigma;

                        //-----------  Tree filling and variables cleaning  ----------//
                    
                        tree_3D->Fill();
                        distances.clear();

                    } else continue;
                } else continue; // quad
            } else continue;
        } // close for pmt
    } // close for cam
    */


    file_root->cd();
    tree_3D->Write();
    file_root->Close();

    deleteNonAlphaDirectories(final_out.c_str());

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Time taken to run the script: " << duration.count() << " seconds" << endl;

    cout << "**Finished**" << endl;
    if(!batch_mode) myapp->Run();

    return 0;
}


/*  *****************************************************************************************************************************************************************************************************************  */

void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E){
    // Function to calculate indices for clusters based on reduced pixels
    // nSc: n_superclusters; npix: total n_pixels
    // sc_redpixID: array of indices for the first reduced pixel in each cluster
    // nSc_red: number of processed superclusters (output)
    // B: begin indices for reduced pixels (output);  E: end indices for reduced pixels (output)

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

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}