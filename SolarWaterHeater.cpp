#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

vector<string> SWH_name = {
	"Solar_intensity", "Factor_conv", "Power_solar_eff",
	"T_water_in_solar", "M_water_solar", "T_water_out_solar"
};

tuple<vector<string>, vector<double>> SWH_Sim(const vector<double>& input, double T_stop) {
	vector<double> para(6);
	para[0] = input[0]; // Solar_intensity (kW/m2)
	para[1] = input[1]; // Factor_conv ( / )
	para[2] = para[0] * para[1]; // Power_solar_eff (kW)
	para[3] = input[2]; // T_water_in (K)

	if (para[3] >= T_stop) {
		para[1] = 0.0;
		para[2] = 0.0;
	}

	para[4] = input[3]; // M_water (kg/s)
	para[5] = para[2] / para[4] / 4.2052 + para[3]; // T_water_out (K), cp_water under 90¡æ

	return { SWH_name,para };
}