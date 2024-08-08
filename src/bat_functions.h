#ifndef BAT_FUNCS_H
#define BAT_FUNCS_H

#include <string>

std::vector<std::string> trim_name(const std::string &full_name, char delimiter);

void coord_change_pmt(double &x, double &y);

void create_bat_input(double run, double event, double trigger, std::vector<std::vector<double>> slice_ints, const std::string& filename = "bat_input.txt");

void run_bat ( const std::string &input_for_bat, const std::string &output_from_bat, const std::string &bat_exe);

void read_bat ( const std::string &output_from_bat, std::vector<std::pair<double,double>>& fitted_points, double &fitted_L, bool verbose);

void calculate_distance(const std::vector<std::pair<double, double>>& points1, const std::vector<std::pair<double, double>>& points2, std::vector<std::pair<double, double>>& distances, bool verbose = false);

#endif