#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_Load.h"

vector<string> ESB_name = {
	"P0", "Load", "P", "M_steam", "T_water_in", "T_steam_out",
	"p_steam_out", "eta_boiler", "h_steam_out"
};

tuple<vector<string>, vector<double>> ESB_Sim(const vector<double>& input) {
	vector<double> para(9);
	para[0] = input[0]; // P0 (kW)
	para[1] = input[1]; // Load (-)
	para[2] = para[0] * para[1]; // P (kW)
	para[3] = input[2]; // M_steam (t/h)
	para[4] = input[3]; // T_water_in (K)
	para[6] = input[4]; // p_steam_out (MPa)
	para[7] = input[5]; // eta_boiler (-)
	para[8] = para[2] * para[7] / (para[3] / 3.6) + 4.1841 * (para[4]-273.15); // h_steam_out (kJ/kg), Cp_water under 20¡æ
	
	para[5] = Cal_T_Steam_P(para[6]); // T_steam_out (K)	

	return { ESB_name,para };
}