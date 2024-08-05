#include <vector>
#include <memory>
#include <stdexcept>
#include <numeric>
#include <iostream>
#include <algorithm>

#include "waveform_analyser.h"

void movingAverageFilter(std::shared_ptr<std::vector<double>>& input_wf, int windowSize) {
    // Ensure the window size is valid
    if (windowSize <= 0) {
        throw std::invalid_argument("Window size must be positive");
    }

    // Get the size of the input_wf vector
    int n = input_wf->size();
    
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

void getTOTs (std::shared_ptr<std::vector<double>> input_wf, double t20_div, double t30_div, 
    int &t20_b, int &t20_e, int &t30_b, int &t30_e, double max, std::vector<int> time) {

    for (int p = 0; p < input_wf->size(); ++p) {
        
        if ( ( (*input_wf)[p] > (max * t20_div) ) && t20_b == 0)                                               t20_b   = time[p];
        if ( ( (*input_wf)[p] < (max * t20_div) ) && t20_b != 0 && t20_e == 0 && ( (time[p] - t20_b) > 10) )   t20_e   = time[p];

        if ( ( (*input_wf)[p] > (max * t30_div) ) && t30_b == 0)                                               t30_b   = time[p];
        if ( ( (*input_wf)[p] < (max * t30_div) ) && t30_b != 0 && t30_e == 0 && ( (time[p] - t30_b) > 10) )   t30_e   = time[p];
    }
    
    if (t20_e == 0) t20_e = 1024;
    if (t30_e == 0) t30_e = 1024;
}

void sliceWaveform_BAT (std::shared_ptr<std::vector<double>> input_wf, 
    std::vector<std::vector<double>> &integrals,
    int nSlices, int TOT20_b, int TOT20_e) {

    double slice_width = (TOT20_e - TOT20_b)/nSlices;
    double peak_int_counter;
    double peak_val;

    integrals.push_back(std::vector<double>{});

    for (int s = 0; s < nSlices; ++s) {

        peak_val = 0;
        peak_int_counter = 0;

        for (int p = TOT20_b+slice_width*(s); p < TOT20_b+slice_width*(s+1); ++p) {
        
            peak_val += (*input_wf)[p];
            peak_int_counter++;
        }

        integrals.back().push_back(peak_val/peak_int_counter);
    }
}

void getSkewness_BraggPeak (std::shared_ptr<std::vector<double>> input_wf, std::vector<int> time_fast_wf,
    int TOT20_b, int TOT20_e, double max_a, int max_t,
    std::vector<double> &skewnesses, std::vector<std::pair<int,double>> &peaks) {

    double TOT20_middle = (TOT20_b + ((TOT20_e-TOT20_b)/2.));

    std::pair<int, double> peak1(0, 0), peak2(0, 0); // (time,amplitude)

    if ( max_t >= TOT20_b && max_t < TOT20_middle ) {       // With this, I look for the "smaller peak" on the opposite half where the max value was found

        for (int p = TOT20_e - 1; p > TOT20_b; --p) {       // Searches for peak starting from the end

            if ( (*input_wf)[p] < (*input_wf)[p+1] && (*input_wf)[p+1] > 0.4*max_a) {       // I know the peak will be in the crown, so I look only for "high amplitude" peaks.

                peak2.first     = time_fast_wf[p+1];
                peak2.second    = (*input_wf)[p+1];
                peak1.first     = max_t;
                peak1.second    = max_a;
                break;
            } 
        }
    
        skewnesses.push_back( (double)( ((-1.)*peak1.second/peak2.second) + 1. ) );
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

        skewnesses.push_back( (double)( ((+1.)*peak2.second/peak1.second) - 1.) );
    }

    peaks.push_back(peak1);
    peaks.push_back(peak2);
}

void getDirectionScore( std::vector<double> skew_ratio, int &dir, double &dir_score, bool verbose) {

    auto    max_skew_it     = std::max_element(skew_ratio.begin(), skew_ratio.end(), [](double a, double b) { return abs(a) < abs(b);});
    double  max_skew        = *max_skew_it;
    int     max_skew_index  = distance(skew_ratio.begin(), max_skew_it);

    for (int i = 0; i < skew_ratio.size(); ++i) {

        if      (i != max_skew_index)   dir_score += skew_ratio[i]/abs(max_skew);  
        else if (i == max_skew_index)   dir_score += (max_skew/abs(max_skew)/2.);   // Need to give some weight to the most clear/obvious waveform. In this case, the weight is 0.5                          // Need to give some weight to the most clear wf. In this case, the weight is 1/2
    }

    if ( verbose == true) {
        std::cout << std::endl;
        std::cout << "--> Bragg Peak Skewness: " << std::endl;                     
        std::cout << "-> Ratios: "; for (auto &val : skew_ratio) std::cout << val << " * "; std::cout<<std::endl;
        std::cout << "-> Abs max ratio: " << abs(max_skew) << std::endl;
        std::cout << "-> Normalized: "; for (auto &val : skew_ratio) std::cout << (double)val/abs(max_skew) << " # "; std::cout<<std::endl;
        std::cout << "** Final Score: " << dir_score*100. << " **" << std::endl;
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

void getAlphaIdentification (std::vector<double> TOT20, std::vector<double> TOT30, bool &pmt_PID_total, bool verbose ) {

    bool pmt_PID1 = false,  pmt_PID2 = false;
    int count_id1 = 0,      count_id2 = 0;

    for (int q = 0; q < TOT20.size(); ++q) {
        
        if ( (TOT20[q] / TOT30[q]) >= 1 && (TOT20[q] / TOT30[q]) <= 2 ) count_id1++;
        
        // if ( TOT20[q] >= 200 ) count_id2++;
        if ( TOT20[q] >= 100 ) count_id2++;
    }

    if ( count_id1 == 4) pmt_PID1 = true;
    if ( count_id2 == 4) pmt_PID2 = true;
    
    if ( pmt_PID1 && pmt_PID2)  pmt_PID_total = true;
    else                        pmt_PID_total = false;

    if (verbose) {

        std::cout << "--> The TOT20/TOT30 ratios were: ";
        for (size_t i = 0; i < TOT20.size(); ++i ) std::cout << static_cast<double>(TOT20[i]/TOT30[i]) << " * ";
        std::cout << std::endl;

        std::cout << "--> The TOT20 lengths were: "; 
        for (const auto &val : TOT20) std::cout << val << " * ";
        std::cout << std::endl;
    }
}