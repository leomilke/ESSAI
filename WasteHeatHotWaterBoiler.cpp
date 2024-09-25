#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

vector<string> WHHWB_name = {
	"Q_in", "eta_heat_recovery", "M_water", "T_water_in", "T_water_out"
};

tuple<vector<string>, vector<double>> WHHWB_Sim(const vector<double> &input) {
	vector<double> para(5);
	para[0] = input[0]; // Q_in (kW)
	para[1] = input[1]; // eta_heat_recovery (-)
	para[2] = input[2]; // M_water (kg/s)
	para[3] = input[3]; // T_water_in (K)

	para[4] = para[0] * para[1] / (4.2156 * para[2]) + para[3]; // T_water_out (K), cp under 100¡æ

	return { WHHWB_name,para };
}