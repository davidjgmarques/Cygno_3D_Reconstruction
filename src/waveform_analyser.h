#ifndef WAVEFORMS_A
#define WAVEFORMS_A

using namespace std;


void movingAverageFilter(std::shared_ptr<std::vector<double>>& input, int windowSize);

void getTOTs (std::shared_ptr<std::vector<double>> input_wf, double t20_div, double t30_div, 
    int &t20_b, int &t20_e, int &t30_b, int &t30_e, double max, std::vector<int> time);

void sliceWaveform_BAT (std::shared_ptr<std::vector<double>> input_wf, 
    std::vector<std::vector<double>> &integrals, int nSlices,
    int TOT20_b, int TOT20_e);

void getQuadrantPMT( std::vector<double>& integrals, int &quadrant_pmt);

void getSkewness_BraggPeak (std::shared_ptr<std::vector<double>> input_wf, std::vector<int> time_fast_wf,
    int TOT20_b, int TOT20_e, double max_a, int max_t,
    std::vector<double> &skewnesses, std::vector<std::pair<int,double>> &peaks);

void getDirectionScore( std::vector<double> skew_ratio, int &dir, double &dir_score, bool verbose);

void getAlphaIdentification (std::vector<double> TOT20, std::vector<double> TOT30,  const int & peaks_ed, bool &pmt_PID_total, bool verbose = true);

void findPeaks(const std::shared_ptr<std::vector<double>>& input_wf, double prominence_percentage, std::vector<std::pair<int, double>>& peaks2, bool verbose = false);


#endif