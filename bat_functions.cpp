#include <vector>
#include <fstream>
#include <iostream>
#include <tuple>
#include <sstream>
#include <algorithm>


#include "bat_functions.h"

// struct bayesd_data {

//     // run_number[0]  event[1]  trigger[2]  peakIndex[3]  L[4]  L_std[5] x[6]  x_std[7]  y[8]  y_std[9]
//     double run, picture, trigger, peakIndex, L, L_std, x, x_std, y, y_std;
// };  

std::vector<std::string> trim_name(const std::string &full_name, char delimiter) {

    std::string trimmed_name;
    std::vector<std::string> tokens;
    std::stringstream check1(full_name);
    std::string intermediate;
    while(getline(check1, intermediate, delimiter)) tokens.push_back(intermediate);

    return tokens;
}

void coord_change_pmt(double &x, double &y) {   //pmt cm to pixels

    double x_new_tmp, y_new_tmp;

    x_new_tmp = x * 1970./33. + 180;
    // y_new_tmp = 2304 - (170 + y*1970./33.);
    y_new_tmp = y * 1970./33. + 370;

    x = x_new_tmp;
    y = y_new_tmp;
}


void create_bat_input(double run, double event, double trigger, std::vector<std::vector<double>> slice_ints, const std::string &filename) {
    
    // We create a file for BAT for eventual debug a posteriori or different studies.
    // In reality, the best would be to simply call a function from bat that does the fit and retrieves specific information directly to the main script.

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

void run_bat ( const std::string &input_for_bat, const std::string &output_from_bat ) {

    std::string bat_executable   = "../BAT_PMTs/./runfit.out";
    std::string bat_input        = " -i " + input_for_bat;
    std::string bat_output       = " -o " + output_from_bat;
    std::string bat_start        = " -s 0";
    std::string bat_end          = " -e 10000";
    std::string bat_method       = " -m association";
    std::string bat_sys_out      = " > bat_system_out.txt";
    std::string bat_command      = bat_executable + bat_input + bat_output + bat_start + bat_end + bat_method + bat_sys_out;
    
    std::cout << std::endl;
    std::cout << "Launching BAT script ..." << std::endl;
    std::cout << bat_command << std::endl;

    int result = system(bat_command.c_str());

    if (result == 0)    std::cout << "-> Script executed successfully." << std::endl;
    else                std::cerr << " !!! Error executing script." << std::endl;
}

void read_bat ( const std::string &output_from_bat, std::vector<std::pair<double,double>>& fitted_points, bool verbose) { 

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

        if (verbose) {
            std::cout <<
                "run: " << run << "; picture: " <<  event << "; trigger: " << trigger   << "; peakIndex: " << peakIndex <<
                "; L: "  << L << "; L_std: " << L_std << "; x: " << x  << "; x_std: " << x_std << "; y: " << y  << "; y_std: " << y_std 
                << std::endl;
        }
    }
    myfile.close();
}

// Function to calculate the X and Y distance between two list of points, making sure that the points are in crescent order of X.
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