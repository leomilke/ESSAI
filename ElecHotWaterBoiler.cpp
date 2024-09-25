#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

vector<string> EHWB_name = {
	"P0", "Load", "P", "M_water", "T_water_in", "T_water_out", "eta_boiler"
};

tuple<vector<string>, vector<double>> EHWB_Sim(const vector<double>& input) {
	vector<double> para(7);
	para[0] = input[0]; // P0 (kW)
	para[1] = input[1]; // Load (-)
	para[2] = para[0] * para[1]; // P (kW)
	para[3] = input[2]; // M_water (t/h)
	para[4] = input[3]; // T_water_in (K)
	para[6] = input[4]; // eta_boiler (-)
	
	para[5] = para[2] * para[6] / (para[3] / 3.6) / 4.2007 + para[4]; // T_water_out (K), Cp_water under 85¡æ	

	return { EHWB_name,para };
}