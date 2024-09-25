#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

vector<string> LiBr_WHW_name = {
	"T_waste_water_in", "M_waste_water", "Q_input", "Q_used",
	"Q_waste_water_out", "T_waste_water_out"
};

tuple<vector<string>, vector<double>> LiBr_WHW_Sim(const vector<double>& input_mix, bool start, vector<double>& x0_pre, string mode) {
	tuple<vector<string>, vector<double>> temp;
	vector<double> para(6);
	vector<double> input;
	input.assign(input_mix.begin() + 2, input_mix.end() - 1);
	double source = *(input_mix.end() - 1);

	if (mode == "Heating") {
		if (source == 0) {
			temp = ASHP_Sim(input, start, x0_pre);
		}
		else if (source == 1) {
			temp = GSHP_Sim(input, mode, start, x0_pre);
		}
		else if (source == 2) {
			temp = WSHP_Sim(input, start, x0_pre);
		}
	}
	else if (mode == "Cooling") {
		if (source == 0) {
			temp = ACC_Sim(input, start, x0_pre);
		}
		else if (source == 1) {
			temp = GSHP_Sim(input, mode, start, x0_pre);
		}
		else if (source == 2) {
			temp = WCC_Sim(input, start, x0_pre);
		}
	}

	vector<string>st = get<0>(temp);
	vector<double>vt = get<1>(temp);

	para[0] = input_mix[0]; // T_waste_water_in (K)
	para[1] = input_mix[1]; // M_waste_water (kg/s)
	para[2] = para[0] * para[1] * 4.2157; // Q_input (kW), Cp_water under 100¡æ
	para[3] = vt[19] / 1.0e3; // Q_used (kW)
	para[4] = para[2] - para[3]; // Q_waste_water_out (kW)
	para[5] = para[4] / para[1] / 4.1850; // T_waste_water_out (K), Cp_water under 60¡æ

	st.insert(st.end(), LiBr_WHW_name.begin(), LiBr_WHW_name.end());
	vt.insert(vt.end(), para.begin(), para.end());

	return{ st,vt };
}