/* -----------------------------------------------------------------------------
 - bat_functions.cpp
 - Small utility helpers for preparing inputs to the external BAT fitter,
 - launching BAT, reading BAT outputs and performing simple coordinate and
 - distance conversions used by the reconstruction/association code.
 - ----------------------------------------------------------------------------- */

#include <vector>
#include <fstream>
#include <iostream>
#include <tuple>
#include <sstream>
#include <algorithm>


#include "bat_functions.h"

/*
 - Split a string by a delimiter and return the tokens as a vector of strings.
 - Parameters:
 -   full_name: input string to be split
 -   delimiter: character delimiting the tokens
 - Returns: vector of string tokens
*/
std::vector<std::string> trim_name(const std::string &full_name, char delimiter) {

    std::string trimmed_name;
    std::vector<std::string> tokens;
    std::stringstream check1(full_name);
    std::string intermediate;
    while(getline(check1, intermediate, delimiter)) tokens.push_back(intermediate);

    return tokens;
}

/*
 - Convert PMT coordinates from cm to detector pixel coordinates.
 - The mapping applied here is an affine transform tuned for this setup
 - and used by downstream plotting/association routines.
*/
void coord_change_pmt(double &x, double &y) {

    double x_new_tmp, y_new_tmp;

    x_new_tmp = x * 1970./33. + 180;
    y_new_tmp = y * 1970./33. + 370;

    x = x_new_tmp;
    y = y_new_tmp;
}

/*
 - Create a BAT input file from PMT slice integrals.
 - Output format (tab-separated): run, event, trigger, slice_index, L1, L2, L3, L4
 - Integrals are converted from ADU to nC using the conversion factor below.
 - Parameters:
 -   run, event, trigger: event identifiers
 -   slice_ints: per-PMT vector of slice integrals (slice_ints.size() = #PMTs)
 -   filename: output filename to write
*/
void create_bat_input(double run, double event, double trigger, std::vector<std::vector<double>> slice_ints, const std::string &filename) {
    
    // The best would be to call a function from bat that does the fit and retrieves specific information directly to the main script.

    /*
     - run: The run number.
     - event: The event number.
     - trigger: The trigger number.
     - peak index: The index indicating the position of the peak in the waveform.
     - L1: The integral of PMT 1.
     - L2: The integral of PMT 2.
     - L3: The integral of PMT 3.
     - L4: The integral of PMT 4.
     - Note: the integral must be given in nC.
     - Conversion: I[nC] = I[ADU] / 4096 * 4/3 * 1/50
    */

    double vtg_to_nC = (1./4096.) * (4./3.) * (1./50.);

    std::ofstream outFile(filename);

    if (outFile.is_open()) {

        for ( int slice_i = 0; slice_i < slice_ints[0].size(); ++slice_i ) {

            outFile << run << "\t" << event << "\t" << trigger << "\t" << slice_i;
            
            for ( int pmt_i = 0; pmt_i < slice_ints.size(); ++pmt_i) {

                outFile << "\t" << slice_ints[pmt_i][slice_i]*vtg_to_nC;
            }
            outFile << "\n";
        }
        outFile.close();
    }
}

/*
 - Run the external BAT executable with the constructed command line.
 - The BAT stdout/stderr are redirected to bat_files/bat_system_out.txt.
 - Parameters:
 -   input_for_bat: path to the input file created by create_bat_input
 -   output_from_bat: path where BAT will write its output
 -   bat_exe: path to the BAT executable
*/
void run_bat ( const std::string &input_for_bat, const std::string &output_from_bat, const std::string &bat_exe) {

    std::string bat_executable   = bat_exe;
    std::string bat_input        = " -i " + input_for_bat;
    std::string bat_output       = " -o " + output_from_bat;
    std::string bat_start        = " -s 0";
    std::string bat_end          = " -e 10000";
    std::string bat_method       = " -m association";
    std::string bat_sys_out      = " > bat_files/bat_system_out.txt";
    std::string bat_command      = bat_executable + bat_input + bat_output + bat_start + bat_end + bat_method + bat_sys_out;
    
    std::cout << std::endl;
    std::cout << "Launching BAT script ..." << std::endl;
    std::cout << bat_command << std::endl;

    int result = system(bat_command.c_str());

    if (result == 0)    std::cout << "-> Script executed successfully." << std::endl;
    else                std::cerr << " !!! Error executing script." << std::endl;
}

/*
 - Read BAT output file and append fitted points to fitted_points.
 - Expects tab-separated columns: run, event, trigger, peakIndex, L, L_std, x, x_std, y, y_std
 - For each parsed line, x and y are transformed to pixel coords via coord_change_pmt
 - fitted_L receives the sum of the L values found in the file
 - Parameters:
 -   output_from_bat: path to BAT output file
 -   fitted_points: output vector (appended)
 -   fitted_L: output parameter containing summed L
 -   verbose: print parsed lines when true
*/
void read_bat ( const std::string &output_from_bat, std::vector<std::pair<double,double>>& fitted_points, double &fitted_L, bool verbose) { 

    double lu = 0;
    double x_corr, y_corr;

    // std::vector<bayesd_data> data;

    std::string line;
    std::ifstream myfile;

    double run, event, trigger, peakIndex, L, L_std, x, x_std, y, y_std;

    myfile.open(output_from_bat);
    std::vector<std::string> name_trim = trim_name(output_from_bat, '/');
    // std::cout << "\nFile opened: " << name_trim.back() << std::endl;

    while ( getline ( myfile, line ) ) {

        std::istringstream iss ( line );	//creates std::string consisting of a line
        std::string token;

        getline (iss, token, '\t'); run         = static_cast<double>(stod(token));         
        getline (iss, token, '\t'); event       = static_cast<double>(stod(token));       
        getline (iss, token, '\t'); trigger     = static_cast<double>(stod(token));     
        getline (iss, token, '\t'); peakIndex   = static_cast<double>(stod(token));   
        getline (iss, token, '\t'); L           = static_cast<double>(stod(token));           
        getline (iss, token, '\t'); L_std       = static_cast<double>(stod(token));       
        getline (iss, token, '\t'); x           = static_cast<double>(stod(token));           
        getline (iss, token, '\t'); x_std       = static_cast<double>(stod(token));       
        getline (iss, token, '\t'); y           = static_cast<double>(stod(token));           
        getline (iss, token, '\t'); y_std       = static_cast<double>(stod(token));       

        // data.push_back( {run, event, trigger, peakIndex, L, L_std, x, x_std, y, y_std} );

        coord_change_pmt(x,y);
        fitted_points.emplace_back(x,y);

        lu += L;

        if (verbose) {
            std::cout <<
                "run: " << run << "; picture: " <<  event << "; trigger: " << trigger   << "; peakIndex: " << peakIndex <<
                "; L: "  << L << "; L_std: " << L_std << "; x: " << x  << "; x_std: " << x_std << "; y: " << y  << "; y_std: " << y_std 
                << std::endl;
        }
    }

    fitted_L = lu;

    myfile.close();
}

/*
 - Compute X/Y distances between two lists of points after sorting both by X.
 - The function requires the two input vectors to have the same length.
 - Parameters:
 -   points1, points2: input lists of (x,y) points
 -   distances: output vector to which (dx,dy) pairs are appended
 -   verbose: if true, print each comparison and resulting distance
*/
void calculate_distance(const std::vector<std::pair<double, double>>& points1, const std::vector<std::pair<double, double>>& points2, std::vector<std::pair<double, double>>& distances, bool verbose) {
    
    if (points1.size() != points2.size()) {
        std::cerr << "Error: The input vectors must have the same size." << std::endl;
        return;
    }

    std::vector<std::pair<double, double>> ordered_points1 = points1;
    std::vector<std::pair<double, double>> ordered_points2 = points2;

    std::sort(ordered_points1.begin(), ordered_points1.end(), [](const std::pair<double, double>& p1, const std::pair<double, double>& p2) {
        return p1.first < p2.first;
    });

    std::sort(ordered_points2.begin(), ordered_points2.end(), [](const std::pair<double, double>& p1, const std::pair<double, double>& p2) {
        return p1.first < p2.first;
    });

    for (int i = 0; i < ordered_points1.size(); ++i) {
        double x_dist = ordered_points2[i].first - ordered_points1[i].first;
        double y_dist = ordered_points2[i].second - ordered_points1[i].second;
        distances.emplace_back(x_dist, y_dist);

        if (verbose) {
            std::cout << "Point 1: (" << ordered_points1[i].first << ", " << ordered_points1[i].second << ")" << std::endl;
            std::cout << "Point 2: (" << ordered_points2[i].first << ", " << ordered_points2[i].second << ")" << std::endl;
            std::cout << "Distance: (" << x_dist << ", " << y_dist << ")" << std::endl;
        }
    }
}