#ifndef Sub_Models
#define Sub_Models
#include <vector>
#include <tuple>
#include <sstream>
#include <numeric>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <algorithm>
#include <functional>
#include <time.h>
#include <thread>  
#include <mutex>
#include <future>
#include <map>
#include <deque>
using namespace std;
//LiBr using hot water

tuple<vector<string>, vector<double>> WSHP_Sim(const vector<double>& input, bool start, vector<double> x0_pre);

tuple<vector<string>, vector<double>> ASHP_Sim(const vector<double>& input, bool start, vector<double> x0_pre);

tuple<vector<string>, vector<double>> WHub_Sim(const vector<double>& T_in, const vector<double>& M_in);

tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> HWSystemSim_single_step(double dt, const vector<string>& heater, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, const vector<double>& T0, vector<int> charge, bool discharge, bool start, vector<vector<double>> x0_pre, vector<vector<double>> FGWasteHeatUse, vector<vector<double>> input_WHHWB);

tuple<vector<string>, vector<double>> WCC_Sim(const vector<double>& input, bool start, vector<double> x0_pre);

tuple<vector<string>, vector<double>> ACC_Sim(const vector<double>& input, bool start, vector<double> x0_pre);

tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> CWSystemSim_single_step(double dt, const vector<string>& chiller, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, const vector<double>& T0, vector<int> charge, bool discharge, bool start, vector<vector<double>> x0_pre);

tuple<vector<string>, vector<double>> GHE_Sim(const vector<double>& input, double M, double T_in, double dh_num, string mode);

tuple<vector<string>, vector<double>> GSHP_Sim(const vector<double>& input_0, string mode, bool start, vector<double> x0_pre);

tuple<vector<string>, vector<double>, vector<double>> STank_Sim(const vector<double>& input, const vector<double>& T0, double dt, string mode);

tuple<vector<double>, vector<tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>>>> HVAC_Sim(double T_total, double dt, const vector<string>& device, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, vector<double> T0, vector<int> charge, bool discharge, int mode, vector<tuple<string, vector<int>>> LB_FG_Source_In, vector<string> LBFG_TP, vector<vector<double>> FGWasteHeatUse, vector<vector<double>> input_WHHWB);

tuple<vector<string>, vector<double>> CBHWB_Sim(const vector<double>& input, bool start, vector<double> T_b_b0);

tuple<vector<string>, vector<double>> CBHWB_Simple_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> GHWB_Sim(const vector<double>& input, bool start, vector<double> T_b_b0);

tuple<vector<string>, vector<double>> GHWB_Simple_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> EHWB_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> CFBHWB_Sim(const vector<double>& input, bool start, vector<double> T_b_b0);

tuple<vector<string>, vector<double>> CFBHWB_Simple_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> ESB_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> SteamHub_Sim(const vector<double>& T_in, const vector<double>& H_in, const vector<double>& M_in, const vector<double>& P_in);

tuple<vector<double>, vector<tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>>>> Steam_Sys_Sim(double T_total, double dt, const vector<string>& device, const vector<vector<double>>& input, map<int, int>& LB_St_Con, int LB_num, vector<tuple<string, vector<int>>> LB_FG_Source_In, vector<string> LBFG_TP, vector<vector<double>> FGWasteHeatUse_S, vector<vector<double>> input_WHSB);

tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>> SteamSysSim_single_step(const vector<string>& heater, const vector<vector<double>>& input, bool start, vector<vector<double>> x0_pre, map<int, int>& LB_St_Con, int LB_num, vector<vector<double>> FGWasteHeatUse_S, vector<vector<double>> input_WHSB);

tuple<vector<string>, vector<double>> CBSB_Sim(const vector<double>& input, bool start, vector<double> T_b_b0);

tuple<vector<string>, vector<double>> CBSB_Simple_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> CFBSB_Sim(const vector<double>& input, bool start, vector<double> T_b_b0);

tuple<vector<string>, vector<double>> CFBSB_Simple_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> GSB_Sim(const vector<double>& input, bool start, vector<double> T_b_b0);

tuple<vector<string>, vector<double>> GSB_Simple_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> LiBr_Steam_Sim(const vector<double>& input, bool start, vector<double>& x0_pre, vector<double>& Steam_in_out, string mode, string source);

tuple<vector<string>, vector<double>> LiBr_DC_Sim(const vector<double>& input_mix, bool start, vector<double>& x0_pre, string mode);

tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> HWSystemSim_single_step_steam(double dt, const vector<string>& heater, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, const vector<double>& T0, vector<int> charge, bool discharge, bool start, vector<vector<double>> x0_pre, vector<vector<double>>& Steam_in_out, vector<string> LB_TP, vector<vector<double>> FGWasteHeatUse, vector<vector<double>> input_WHHWB);

tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> CWSystemSim_single_step_steam(double dt, const vector<string>& chiller, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, const vector<double>& T0, vector<int> charge, bool discharge, bool start, vector<vector<double>> x0_pre, vector<vector<double>>& Steam_in_out, vector<string> LB_TP);

tuple<vector<double>, vector<tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>>>, vector<tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>>>> HVAC_Steam_Sim(double T_total, double dt, const vector<string>& HVAC_device, const vector<string>& Steam_device, const vector<vector<double>>& input_HVAC, const vector<vector<double>>& input_Steam, bool with_storage, const vector<double>& input_ST, vector<double> T0, vector<int> charge, bool discharge, int mode, vector<string>& LB_Type, vector<vector<int>>& LB_St_Con, vector<tuple<string, vector<int>>> LB_FG_Source_In, vector<string> LBFG_TP, vector<vector<double>> FGWasteHeatUse_W, vector<vector<double>> input_WHHWB, vector<vector<double>> FGWasteHeatUse_S, vector<vector<double>> input_WHSB);

tuple<vector<string>, vector<double>> FGHub_Sim(const vector<double>& Q_in, const vector<double>& T_in);

tuple<vector<string>, vector<double>> LiBr_FG_Sim(const vector<double>& input, bool start, vector<double>& x0_pre, string mode, string source, vector<double> FG_input);

tuple<vector<string>, vector<double>> CTowerOpen_Sim(const vector<double>& input, bool start, bool iter_start);

tuple<vector<string>, vector<double>> CTowerOpen_Simple_Sim(const vector<double>& input, bool start, bool iter_start);

tuple<vector<string>, vector<double>> CTowerClosed_Simple_Sim(const vector<double>& input, bool start, bool iter_start);

tuple<vector<string>, vector<double>> CTowerClosed_Sim(const vector<double>& input, bool start, bool iter_start);

tuple<vector<string>, vector<double>> WHHWB_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> WHSB_Sim(const vector<double>& input);

tuple<vector<string>, vector<double>> LiBr_WHW_Sim(const vector<double>& input_mix, bool start, vector<double>& x0_pre, string mode);

tuple<vector<string>, vector<double>> SWH_Sim(const vector<double>& input, double T_stop);

tuple<vector<string>, vector<double>> LiBr_Solar_Sim(const vector<double>& input_mix, bool start, vector<double>& x0_pre, string mode, double dt);

#endif