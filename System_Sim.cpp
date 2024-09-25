#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_Load.h"

// In Gas Boiler, ¦Á and p_steam_out should only vary with Load.

deque<double> Load_pre;
deque<vector<double>> K_combust_boiler, K_CoolingTower, ST_T_inside;
deque<string> LB_FG_TYPE;
deque<int> LB_FG_INDEX;
vector<tuple<string, vector<int>>> LB_FG_Source;
vector<vector<double>> FGInput_Steam, FGInput_Water;

void FlueGasLineInit(vector<string> Steam_Dev, vector<string> HVAC_Dev, vector<tuple<string, vector<int>>> LB_FG_Source_In, vector<string> TP) {
	LB_FG_Source = LB_FG_Source_In;
	vector<vector<double>> Input_Steam(Steam_Dev.size(), vector<double>(2, 0.0));
	FGInput_Steam = Input_Steam;
	vector<vector<double>> Input_Water(HVAC_Dev.size(), vector<double>(2, 0.0));
	FGInput_Water = Input_Water;
	int num_LBFG = count(HVAC_Dev.begin(), HVAC_Dev.end(), "LiBr_FG");
	for (int i = 0; i < num_LBFG; ++i) {
		LB_FG_INDEX.push_back(i);
		LB_FG_TYPE.push_back(TP[i]);
	}
}

tuple<vector<double>, vector<tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>>>, vector<tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>>>> HVAC_Steam_Sim(double T_total, double dt, const vector<string>& HVAC_device, const vector<string>& Steam_device, const vector<vector<double>>& input_HVAC, const vector<vector<double>>& input_Steam, bool with_storage, const vector<double>& input_ST, vector<double> T0, vector<int> charge, bool discharge, int mode, vector<string>& LB_Type, vector<vector<int>>& LB_St_Con, vector<tuple<string, vector<int>>> LB_FG_Source_In, vector<string> LBFG_TP, vector<vector<double>> FGWasteHeatUse_W, vector<vector<double>> input_WHHWB, vector<vector<double>> FGWasteHeatUse_S, vector<vector<double>> input_WHSB) {
	Read_Steam_Dat();
	FlueGasLineInit(Steam_device, HVAC_device, LB_FG_Source_In, LBFG_TP);
	int time_step_num = T_total / dt;
	int num_HVAC_device = HVAC_device.size();
	int num_Steam_device = Steam_device.size();

	vector<tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>>> R(time_step_num);
	vector<tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>>> Rs(time_step_num);
	vector<double> Time(time_step_num);
	tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> temp;
	tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>> temp_Steam;
	vector<vector<double>> x0_pre_HVAC(num_HVAC_device);
	vector<vector<double>> x0_pre_Steam(num_Steam_device);
	vector<vector<double>> St2H (num_HVAC_device);

	map<int, int> LB_St_map;

	vector<int> LiBr_pos;
	vector<string> LiBr_Type(num_HVAC_device);
	for (int p = 0; p < num_HVAC_device; ++p) {
		if (HVAC_device[p] == "LiBr_S") {
			LiBr_pos.push_back(p);
		}
	}
	int LiBr_num = LiBr_pos.size();
	for (int p = 0; p < LiBr_num; ++p) {
		LiBr_Type[LiBr_pos[p]] = LB_Type[p];
	}

	for (int p = 0; p < LiBr_num; ++p) {
		for (int pp : LB_St_Con[p]) { LB_St_map.insert({ pp,p }); }
	}
	
	if (mode == 0) { // Heating
		for (int i = 0; i < time_step_num; ++i) {
			temp_Steam = SteamSysSim_single_step(Steam_device, input_Steam, !i, x0_pre_Steam, LB_St_map, LiBr_num, FGWasteHeatUse_S, input_WHSB);
			for (int j = 0; j < LiBr_num; ++j) { St2H[LiBr_pos[j]] = get<3>(temp_Steam)[j]; }
			temp = HWSystemSim_single_step_steam(dt, HVAC_device, input_HVAC, with_storage, input_ST, T0, charge, discharge, !i, x0_pre_HVAC, St2H, LiBr_Type, FGWasteHeatUse_W, input_WHHWB);
			R[i] = temp;
			Rs[i] = temp_Steam;
			T0 = get<2>(temp);
			x0_pre_HVAC = get<3>(temp);
			x0_pre_Steam = get<2>(temp_Steam);
			Time[i] = ((i + 1.0) * dt);
		}
	}
	else if (mode == 1) { // Cooling
		for (int i = 0; i < time_step_num; ++i) {
			temp_Steam = SteamSysSim_single_step(Steam_device, input_Steam, !i, x0_pre_Steam, LB_St_map, LiBr_num, FGWasteHeatUse_S, input_WHSB);
			for (int j = 0; j < LiBr_num; ++j) { St2H[LiBr_pos[j]] = get<3>(temp_Steam)[j]; }
			temp = CWSystemSim_single_step_steam(dt, HVAC_device, input_HVAC, with_storage, input_ST, T0, charge, discharge, !i, x0_pre_HVAC, St2H, LiBr_Type);
			R[i] = temp;
			Rs[i] = temp_Steam;
			T0 = get<2>(temp);
			x0_pre_HVAC = get<3>(temp);
			x0_pre_Steam = get<2>(temp_Steam);
			Time[i] = ((i + 1.0) * dt);
		}
	}

	return{ Time,R,Rs };
}

tuple<vector<double>, vector<tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>>>> HVAC_Sim(double T_total, double dt, const vector<string>& device, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, vector<double> T0, vector<int> charge, bool discharge, int mode, vector<tuple<string, vector<int>>> LB_FG_Source_In, vector<string> LBFG_TP, vector<vector<double>> FGWasteHeatUse_W, vector<vector<double>> input_WHHWB) {
	FlueGasLineInit({}, device, LB_FG_Source_In, LBFG_TP);
	int time_step_num = T_total / dt;
	vector<tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>>> R(time_step_num);
	vector<double> Time(time_step_num);
	tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> temp;
	vector<vector<double>> x0_pre(device.size());

	if (mode == 0) { // Heating
		for (int i = 0; i < time_step_num; ++i) {
			temp = HWSystemSim_single_step(dt, device, input, with_storage, input_ST, T0, charge, discharge, !i, x0_pre, FGWasteHeatUse_W, input_WHHWB);
			R[i] = temp;
			T0 = get<2>(temp);
			x0_pre = get<3>(temp);
			Time[i] = ((i + 1.0) * dt);
		}
	}
	else if (mode == 1) { // Cooling
		for (int i = 0; i < time_step_num; ++i) {
			temp = CWSystemSim_single_step(dt, device, input, with_storage, input_ST, T0, charge, discharge, !i, x0_pre);
			R[i] = temp;
			T0 = get<2>(temp);
			x0_pre = get<3>(temp);
			Time[i] = ((i + 1.0) * dt);
		}
	}

	return{ Time,R };
}

tuple<vector<double>, vector<tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>>>> Steam_Sys_Sim(double T_total, double dt, const vector<string>& device, const vector<vector<double>>& input, map<int, int>& LB_St_Con, int LB_num, vector<tuple<string, vector<int>>> LB_FG_Source_In, vector<string> LBFG_TP, vector<vector<double>> FGWasteHeatUse_S, vector<vector<double>> input_WHSB) {
	Read_Steam_Dat();
	FlueGasLineInit(device, {}, LB_FG_Source_In, LBFG_TP);
	int time_step_num = T_total / dt;
	vector<tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>>> R(time_step_num);
	vector<double> Time(time_step_num);
	tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>> temp;
	vector<vector<double>> x0_pre(device.size());

	for (int i = 0; i < time_step_num; ++i) {
		temp = SteamSysSim_single_step(device, input, !i, x0_pre, LB_St_Con, LB_num, FGWasteHeatUse_S, input_WHSB);
		R[i] = temp;
		x0_pre = get<2>(temp);
		Time[i] = ((i + 1.0) * dt); 
	}

	return{ Time,R };
}