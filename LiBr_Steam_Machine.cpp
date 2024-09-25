#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

vector<string> LiBr_Steam_name = {
	"H_steam_in", "M_steam", "Q_steam_in", "Q_steam_used",
	"Q_con_water_out", "H_con_water_out", "T_con_water"
};
// Connect to WHSB?
tuple<vector<string>, vector<double>> LiBr_Steam_Sim(const vector<double>& input, bool start, vector<double>& x0_pre, vector<double>& Steam_in_out, string mode, string source) {
	tuple<vector<string>, vector<double>> temp;
	vector<double> para(7);

	if (mode == "Heating") {
		if (source == "Air") {
			temp = ASHP_Sim(input, start, x0_pre);
		}
		else if (source == "Ground") {
			temp = GSHP_Sim(input, mode, start, x0_pre);
		}
		else if (source == "Water") {
			temp = WSHP_Sim(input, start, x0_pre);
		}
	}
	else if (mode == "Cooling") {
		if (source == "Air") {
			temp = ACC_Sim(input, start, x0_pre);
		}
		else if (source == "Ground") {
			temp = GSHP_Sim(input, mode, start, x0_pre);
		}
		else if (source == "Water") {
			temp = WCC_Sim(input, start, x0_pre);
		}
	}
	
	vector<string>st = get<0>(temp);
	vector<double>vt = get<1>(temp);

	para[0] = Steam_in_out[0]; // H_steam_in (kJ/kg)
	para[1] = Steam_in_out[1]; // M_steam (t/h)
	para[2] = para[0] * (para[1] / 3.6); // Q_steam_in (kW)
	para[3] = vt[19] / 1.0e3; // Q_steam_used (kW)
	para[4] = para[2] - para[3]; // Q_con_water_out (kW)
	para[5] = para[4] / (para[1] / 3.6); // H_con_water_out (kJ/kg)
	para[6] = para[5] / 4.2007 + 273.15; // T_con_water (K), Cp_water under 85¡æ

	st.insert(st.end(), LiBr_Steam_name.begin(), LiBr_Steam_name.end());
	vt.insert(vt.end(), para.begin(), para.end());

	return{ st,vt };
}