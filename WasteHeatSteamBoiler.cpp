#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_Load.h"

vector<string> WHSB_name = {
	"Q_in", "eta_heat_recovery", "M_water/steam", "T_water_in", "h_steam_out", "P_steam_out", "T_steam_out"
};

tuple<vector<string>, vector<double>> WHSB_Sim(const vector<double>& input) {
	vector<double> para(7);
	para[0] = input[0]; // Q_in (kW)
	para[1] = input[1]; // eta_heat_recovery (-)
	para[2] = input[2]; // M_water/steam (kg/s)
	para[3] = input[3]; // T_water_in (K)

	para[4] = para[0] * para[1] / para[2] + 4.1841 * para[3]; // h_steam_out (kJ/kg), Cp_water under 20¡æ
	para[5] = input[4]; // p_steam_out (MPa)
	para[6] = Cal_T_Steam_P(para[5]); // T_steam_out (K)

	return { WHSB_name,para };
}