#ifndef BAT_FUNCS_H
#define BAT_FUNCS_H

// struct bayesd_data;

std::vector<std::string> trim_name(const std::string &full_name, char delimiter);
void coord_change_pmt(double &x, double &y);

void create_bat_input(double run, double event, double trigger, std::vector<std::vector<double>> slice_ints, const std::string& filename = "bat_input.txt");
void run_bat ( const std::string &input_for_bat, const std::string &output_from_bat );
void read_bat ( const std::string &output_from_bat, std::vector<std::pair<double,double>>& fitted_points, bool verbose = false);

#endif