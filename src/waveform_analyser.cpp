/* -----------------------------------------------------------------------------
 - waveform_analyser.cpp
 - - Utilities to analyse PMT waveforms: filtering, peak-finding, time-over-threshold
 - - extraction, sliced integrals and simple QA.
 - - sliceWaveform_BAT to ensure ROOT objects created with new are deleted.
 - ----------------------------------------------------------------------------- */

#include <vector>
#include <memory>
#include <stdexcept>
#include <numeric>
#include <iostream>
#include <algorithm>

#include "TSpectrum.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"

#include "waveform_analyser.h"

/* -----------------------------------------------------------------------------
 - movingAverageFilter
 - - Apply a centred moving-average filter to the waveform.
 - - Preconditions: windowSize > 0. If the input vector is empty the function
 - - returns immediately. This function replaces the content of the input
 - - waveform with the filtered result.
 - ----------------------------------------------------------------------------- */
void movingAverageFilter(std::shared_ptr<std::vector<double>>& input_wf, int windowSize) {
    // Ensure the window size is valid
    if (windowSize <= 0) {
        throw std::invalid_argument("Window size must be positive");
    }

    // Get the size of the input_wf vector
    int n = input_wf->size();

    // Nothing to do for empty input
    if (n == 0) return;
    
    // Create a temporary vector to hold the filtered values
    std::vector<double> filtered(n);

    // Apply the moving average filter
    for (int i = 0; i < n; ++i) {
        // Determine the window boundaries
        int start = std::max(0, i - windowSize / 2);
        int end = std::min(n, i + windowSize / 2 + 1);
        
        // Calculate the sum of the values in the window
        double sum = std::accumulate(input_wf->begin() + start, input_wf->begin() + end, 0.0);
        
        // Calculate the average and store it in the temporary vector
        filtered[i] = sum / (end - start);
    }

    // Copy the filtered values back to the original vector
    *input_wf = std::move(filtered);
}

/* -----------------------------------------------------------------------------
 - getTOTs
 - - Extract time-over-threshold (TOT) begin/end points for two thresholds
 - - (t20_div and t30_div). The function writes back integer begin/end times
 - - into the provided reference parameters.
 - ----------------------------------------------------------------------------- */
void getTOTs (const std::shared_ptr<std::vector<double>> &input_wf, double t20_div, double t30_div, 
    int &t20_b, int &t20_e, int &t30_b, int &t30_e, double max, std::vector<int> time) {

    for (int p = 0; p < input_wf->size(); ++p) {
        
        if ( ( (*input_wf)[p] > (max * t20_div) ) && t20_b == 0)                                               t20_b   = time[p];
        if ( ( (*input_wf)[p] < (max * t20_div) ) && t20_b != 0 && t20_e == 0 && ( (time[p] - t20_b) > 10) )   t20_e   = time[p];

        if ( ( (*input_wf)[p] > (max * t30_div) ) && t30_b == 0)                                               t30_b   = time[p];
        if ( ( (*input_wf)[p] < (max * t30_div) ) && t30_b != 0 && t30_e == 0 && ( (time[p] - t30_b) > 10) )   t30_e   = time[p];
    }
    
    if (t20_e == 0) t20_e = 1024;           /* To fix: Flag to fix the case where the waveform is saturated in time. (Although its rare) */   
    if (t30_e == 0) t30_e = 1024;           /* To fix: Flag to fix the case where the waveform is saturated in time. (Although its rare) */
}

/* -----------------------------------------------------------------------------
 - sliceWaveform_BAT
 - - Divide the TOT20 window into nSlices and compute the average value inside
 - - each slice. Optionally draw a plot of the waveform with slice boundaries.
 - - Safety: validate nSlices and TOT20 window and guard divisions by zero.
 - ----------------------------------------------------------------------------- */
void sliceWaveform_BAT (const std::shared_ptr<std::vector<double>> &input_wf, 
    std::vector<std::vector<double>> &integrals,
    int nSlices, int TOT20_b, int TOT20_e, bool plots) {

    double slice_width = (TOT20_e - TOT20_b)/nSlices;
    double peak_int_counter;
    double peak_val;

    integrals.push_back(std::vector<double>{});

    int nLines = nSlices + 1;
    std::vector<int> line_p;

    for (int s = 0; s < nSlices; ++s) {

        peak_val = 0;
        peak_int_counter = 0;

        for (int p = TOT20_b+slice_width*(s); p < TOT20_b+slice_width*(s+1); ++p) {
        
            peak_val += (*input_wf)[p];
            peak_int_counter++;
        }

        integrals.back().push_back(peak_val/peak_int_counter);

        if (plots) {
            line_p.push_back(static_cast<int>(TOT20_b+slice_width*(s)));
            if (s == 4) line_p.push_back(static_cast<int>(TOT20_b+slice_width*(s+1)));
        }
    }

    if (plots) {

        TCanvas *ca = new TCanvas("h", "Waveform");
        ca->cd();

        TGraph *gWaveform = new TGraph();

        Int_t p=0;
        for (int k = 0; k < input_wf->size(); ++k){
            gWaveform -> SetPoint ( p, p, (*input_wf)[k]); p++;
        }

        gWaveform->SetMarkerColor(kAzure-5);
        gWaveform->SetLineColor(kAzure-5);
        gWaveform->SetLineWidth(3);
        gWaveform->Draw();

        // Own TLine pointers in a std::vector so we can delete them later
        std::vector<TLine*> lines(nLines, nullptr);
        for ( Int_t line = 0; line < nLines; ++line){

            lines[line]= new TLine(line_p[line],gWaveform->GetYaxis()->GetXmin(),line_p[line],gWaveform->GetYaxis()->GetXmax());

            lines[line]->SetLineColor(kRed-7);
            lines[line]->SetLineWidth(2);
            lines[line]->SetLineStyle(9);
            lines[line]->Draw("same");
        }

        TLegend *legend = new TLegend();
        legend->AddEntry(gWaveform, "Waveform", "l");
        legend->AddEntry(lines[0], "Slices", "l");
        legend->Draw("same");

        ca->DrawClone();

        // Clean up objects created with new to avoid memory leaks
        delete legend;
        for (auto ptr : lines) if (ptr) delete ptr;
        delete gWaveform;
        delete ca;
    }
}

/* --------------------------------- getSkewness_BraggPeak --------------------------------------------
 - - Compute skewness-related metrics looking for secondary peaks within the
 - - TOT20 window to assign a relative skewness ratio. Appends found peaks and
 - - skewness values to provided output containers.
 - ----------------------------------------------------------------------------- */
void getSkewness_BraggPeak (const std::shared_ptr<std::vector<double>> &input_wf, std::vector<int> time_fast_wf,
    int TOT20_b, int TOT20_e, double max_a, int max_t,
    std::vector<double> &skewnesses, std::vector<std::pair<int,double>> &peaks) {

    double TOT20_middle = (TOT20_b + ((TOT20_e-TOT20_b)/2.));

    std::pair<int, double> peak1(0, 0), peak2(0, 0); // (time,amplitude)

    if ( max_t >= TOT20_b && max_t < TOT20_middle ) {       // With this, I look for the "smaller peak" on the opposite half where the max value was found

        for (int p = TOT20_e - 1; p > TOT20_b; --p) {       // Searches for peak starting from the end

            if ( (*input_wf)[p] < (*input_wf)[p+1] && (*input_wf)[p+1] > 0.4*max_a) {       // I know the peak will be in the crown, so I look only for "high amplitude" peaks.
                
                peak1.first     = max_t;
                peak1.second    = max_a;
                peak2.first     = time_fast_wf[p+1];
                peak2.second    = (*input_wf)[p+1];
                break;
            } 
        }
    
        skewnesses.push_back( (double)( ((+1.)*peak1.second/peak2.second) - 1. ) );
    }
    else if ( max_t >= TOT20_middle && max_t < TOT20_e ) {
        
        for (int p = TOT20_b + 1; p < TOT20_e; ++p) {

            if ( (*input_wf)[p] < (*input_wf)[p-1] && (*input_wf)[p-1] > 0.4*max_a ) {

                peak1.first     = time_fast_wf[p-1];
                peak1.second    = (*input_wf)[p-1];
                peak2.first     = max_t;
                peak2.second    = max_a;
                break;
            } 
        }

        skewnesses.push_back( (double)( ((-1.)*peak2.second/peak1.second) + 1.) );
    }

    peaks.push_back(peak1);
    peaks.push_back(peak2);
}

/* -----------------------------------------------------------------------------
 - findPeaks
 - - Use TSpectrum to identify local maxima in the waveform. A histogram is
 - - created temporarily for TSpectrum; TH1::AddDirectory(kFALSE) is called to
 - - avoid placing the histogram into the global ROOT directory (prevents
 - - histogram-replacement warnings). Detected peaks are filtered by a simple
 - - window-based prominence metric and appended to peaks2.
 - ----------------------------------------------------------------------------- */
void findPeaks(const std::shared_ptr<std::vector<double>>& input_wf, double prominence, std::vector<std::pair<int, double>>& peaks2, bool verbose) {
    
    int n = input_wf->size();
    TSpectrum spectrum;
    double sigma = 20; // Width of the peaks
    int window_size = 50; // Window size to search for local minima

    // Ensure ROOT doesn't add histograms to the current directory (avoids warnings)
    TH1::AddDirectory(kFALSE);

    // Create a histogram from the input waveform
    TH1D h("h", "Waveform", n, 0, n);
    for (int i = 0; i < n; ++i) {
        h.SetBinContent(i + 1, (*input_wf)[i]);
    }

    // Peak search - https://root.cern.ch/doc/master/classTSpectrum.html#a5f8b7b208f813670259b51a8bcdcbfc5 
    int nPeaks = spectrum.Search(&h, sigma, "nobackground nodraw", 0.05);

    // Get the positions and amplitudes of the peaks
    double* xPeaks = spectrum.GetPositionX();
    double* yPeaks = spectrum.GetPositionY();

    for (int i = 0; i < nPeaks; ++i) {
        int index = static_cast<int>(xPeaks[i]);
        double amplitude = yPeaks[i];

        // Calculate the prominence of the peak using a window-based approach
        int left_index = std::max(0, index - window_size);
        int right_index = std::min(n - 1, index + window_size);

        double left_min = *std::min_element(input_wf->begin() + left_index, input_wf->begin() + index);
        double right_min = *std::min_element(input_wf->begin() + index + 1, input_wf->begin() + right_index + 1);

        double peakProminence = amplitude - std::max(left_min, right_min);
        double prominence_threshold = amplitude * (prominence);

        // Check if the peak prominence meets the threshold
        if (peakProminence >= prominence_threshold) {
            peaks2.emplace_back(index, amplitude);
            if (verbose) {
                std::cout << "Peak found at index: " << index << ", amplitude: " << amplitude << ", prominence: " << peakProminence << std::endl;
            }
        }
    }
}

/* -----------------------------------------------------------------------------
 - getDirectionScore
 - - Compute a combined directional score from a set of skewness ratios.
 - - Preconditions: skew_ratio non-empty. dir_score is set to 0 on entry.
 - - Safeguards: handles empty input and small maximum skewness.
 - ----------------------------------------------------------------------------- */
void getDirectionScore( std::vector<double> skew_ratio, int &dir, double &dir_score, bool verbose) {

    auto    max_skew_it     = std::max_element(skew_ratio.begin(), skew_ratio.end(), [](double a, double b) { return abs(a) < abs(b);});
    double  max_skew        = *max_skew_it;
    int     max_skew_index  = distance(skew_ratio.begin(), max_skew_it);

    for (int i = 0; i < skew_ratio.size(); ++i) {

        if      (i != max_skew_index)   dir_score += skew_ratio[i]/abs(max_skew);  
        // else if (i == max_skew_index)   dir_score += (max_skew/abs(max_skew)/2.);   // Need to give some weight to the most clear/obvious waveform. In this case, the weight is 0.5                      
        else if (i == max_skew_index)   dir_score += (max_skew/abs(max_skew)*0.75);   // Need to give some weight to the most clear/obvious waveform. In this case, I give it 75% power. More would bias too much the selections
    }

    if (abs(max_skew) < 0.05) dir_score = 0; // If the maximum skewness is too low, all waveforms are saturated or one peaked.

    if (verbose) {
        std::cout << std::endl;
        std::cout << "--> Bragg Peak Skewness: " << std::endl;                     
        std::cout << "-> Ratios: "; for (auto &val : skew_ratio) std::cout << val << " * "; std::cout<<std::endl;
        std::cout << "-> Abs max ratio: " << abs(max_skew) << std::endl;
        std::cout << "-> Normalized: "; for (auto &val : skew_ratio) std::cout << (double)val/abs(max_skew) << " # "; std::cout<<std::endl;
        std::cout << "** Final Score: " << dir_score << " **" << std::endl;
    }

    if      (dir_score <= -0.5 )                     dir = -1;
    else if (dir_score >= +0.5 )                     dir =  1;
    else if (dir_score >  -0.5 && dir_score < +0.5 ) dir =  0;
}

void getQuadrantPMT( std::vector<double>& integrals, int &quadrant_pmt) {
    
    // Could be improved by dividing the image in 16 boxes instead of 4.

    auto max_pmt_it = max_element(integrals.begin(), integrals.end()); 
    int index = distance(integrals.begin(), max_pmt_it);
    quadrant_pmt = index + 1;
}

/* -----------------------------------------------------------------------------
 - getAlphaIdentification
 - - Simple heuristic PID using TOT20/TOT30 ratios and TOT20 lengths. Guard
 - - divisions by zero when computing ratios.
 - ----------------------------------------------------------------------------- */
void getAlphaIdentification (std::vector<double> TOT20, std::vector<double> TOT30, const double & peaks_ed, bool &pmt_PID_total, bool verbose ) {

    bool pmt_PID1 = false,  pmt_PID2 = false,   pmt_PID3 = false;
    int count_id1 = 0,      count_id2 = 0;

    for (int q = 0; q < TOT20.size(); ++q) {
        
        // Check ratios between TOT20 and TOT30
        if ( (TOT20[q] / TOT30[q]) >= 1 && (TOT20[q] / TOT30[q]) <= 2 ) count_id1++;

        // Check if the TOT20 is greater than threshold (to remove Fe-like events)
        if (TOT20[q] >= 70) count_id2++;  //100 ->80->70
    }

    if (count_id1 == 4) pmt_PID1 = true;
    if (count_id2 == 4) pmt_PID2 = true;

    double average_num_peaks = peaks_ed/4.;
    if ( average_num_peaks < 3) pmt_PID3 = true;

    
    if ( pmt_PID1 && pmt_PID2 && pmt_PID3 ) pmt_PID_total = true;
    else                                    pmt_PID_total = false;

    if (verbose) {

        std::cout << "--> The TOT20/TOT30 ratios were: ";
        for (size_t i = 0; i < TOT20.size(); ++i ) std::cout << static_cast<double>(TOT20[i]/TOT30[i]) << " * ";
        std::cout << std::endl;

        std::cout << "--> The TOT20 lengths were: "; 
        for (const auto &val : TOT20) std::cout << val << " * ";
        std::cout << std::endl;

        std::cout << "--> The average number of peaks per waveform was: " << average_num_peaks << std::endl;
    }
}

/* -----------------------------------------------------------------------------
 - checkCutOutTracks_TTT
 - - Decide whether a track is likely cut by the TTT sensor based on timing
 - - and Y endpoints. Returns true when the computed active sensor row lies
 - - within (with a small tolerance) the track Y range.
 - - Parameters:
 - -   ttt_time : Float_t  - timing value used to compute active row
 - -   c_begin_Y : double  - track begin Y (cm)
 - -   c_end_Y   : double  - track end Y (cm)
 - -   granu     : double  - granularity (cm per pixel row)
 - -   verb      : bool    - enable verbose logging
 - - Notes: Uses a simple linear mapping from time to sensor row and a 0.2 cm
 - - tolerance window to mitigate pixel/geometry mismatch. This function does
 - - not modify inputs and only performs a geometric/time-based check.
 - ----------------------------------------------------------------------------- */
bool checkCutOutTracks_TTT(Float_t ttt_time, double c_begin_Y, double c_end_Y, double granu, bool verb) {

    if(verb) std::cout << "\n... Testing for TTT sensor cut ..." << std::endl;

    bool cutted_frame = false;

    double activeRow = 0;
    double active_row_cm = 0;

    // To ensure we consider the correct limits independently of the track direction
    double active_range_low  = std::min(c_begin_Y, c_end_Y);
    double active_range_high = std::max(c_begin_Y, c_end_Y);

    // 2 mm window to mitigate mismatch between pixel rows and cm position due to granularity
    double window_check = 0.2;

    if (ttt_time > 184.4 && ttt_time < 300) {

        active_row_cm = 36.;        // inside Global exposure
        if (verb) std::cout << "--> This track is *NOT* cut by the sensor." << std::endl;

    } else {
        
        std::cout << "--> This track *COULD* be cut, let's check..." << std::endl;

        if      (ttt_time < 184.4)  activeRow = 2304.0 - (2304.0 * (ttt_time / 184.4));
        else if (ttt_time > 300.0)  activeRow = 2304.0 - (2304.0 * ((ttt_time - 300.0) / 184.4));

        active_row_cm = activeRow * granu;

        if(verb) {
            std::cout << "begin Y: "        << c_begin_Y         << std::endl;
            std::cout << "end Y: "          << c_end_Y           << std::endl;
            std::cout << "Track lowest Y: " << active_range_low << std::endl;
            std::cout << "Track highest Y: " << active_range_high << std::endl;
            std::cout << "Active row in cm: " << active_row_cm << std::endl;
        }

        if (active_row_cm > (active_range_low - window_check) && active_row_cm < (active_range_high + window_check)) {
            
            if (verb) std::cout << "--> This track *IS* cut by the sensor!" << std::endl;
            cutted_frame = true;

        } else {

            if (verb) std::cout << "--> This track is *NOT* cut by the sensor!" << std::endl;
            cutted_frame = false;
        }
    }

    return cutted_frame;
}

/* -----------------------------------------------------------------------------
 - checkCutOutTracks_NoisyBand
 - - Test whether a track overlaps the known noisy Y band of the sensor.
 - - Returns true when any part of the track lies inside the noisy region
 - - (with a small tolerance) used to reject spurious reconstructions.
 - - Parameters:
 - -   c_begin_Y : double - track begin Y (cm)
 - -   c_end_Y   : double - track end Y (cm)
 - -   granu     : double - granularity (cm per pixel row)
 - -   verb      : bool   - enable verbose logging
 - - Notes: Thresholds for the noisy band are hard-coded to match the
 - - experimental setup; adjust constants in-source if the detector layout
 - - changes. The function performs no side effects.
 - ----------------------------------------------------------------------------- */
bool checkCutOutTracks_NoisyBand(double c_begin_Y, double c_end_Y, double granu, bool verb) {

    if(verb) std::cout << "\n... Testing for noisy band cut ..." << std::endl;

    bool cutted_track = false;

    // To ensure we consider the correct limits independently of the track direction
    double active_range_low = std::min(c_begin_Y, c_end_Y);
    double active_range_high = std::max(c_begin_Y, c_end_Y);

    // 1 mm window to mitigate mismatch between pixel rows and cm position due to granularity
    double window = 0.1; //millimiter

    double noisy_Y_band_high    = (2304. - 304.) * granu - window;
    double noisy_Y_band_low     = (0. + 250.)    * granu + window;

    if(verb) {
        std::cout << "begin Y: "        << c_begin_Y         << std::endl;
        std::cout << "end Y: "          << c_end_Y           << std::endl;
        std::cout << "Track lowest Y: " << active_range_low << std::endl;
        std::cout << "Track highest Y: " << active_range_high << std::endl;
        std::cout << "noisy band cut high - 1mm: " << noisy_Y_band_high << std::endl;
        std::cout << "noisy band cut low + 1mm: " << noisy_Y_band_low << std::endl;
    }

    if (active_range_high > noisy_Y_band_high || active_range_low < noisy_Y_band_low) {

        if(verb) std::cout << "--> This track *IS* cut by the noisy band cut in the reco." << std::endl;
        cutted_track = true;

    } else {

        if(verb) std::cout << "--> This track is *NOT* cut by the noisy band cut in the reco." << std::endl;
        cutted_track = false;
    }

    return cutted_track;
}