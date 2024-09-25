#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_Load.h"

vector<string> LiBr_FG_name = {
	"M_fluegas", "h_fluegas_in", "T_fluegas_in", "Q_fluegas_in", "Q_used",
	"Q_fluegas_out", "h_fluegas_out", "T_fluegas_out"
};

tuple<vector<string>, vector<double>> LiBr_FG_Sim(const vector<double>& input, bool start, vector<double>& x0_pre, string mode, string source, vector<double> FG_input) {
	tuple<vector<string>, vector<double>> temp;
	vector<double> para(8);

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

	para[0] = FG_input[0]; // M_fluegas (kg/s)
	para[1] = FG_input[1]; // h_fluegas_in (kJ/kg)
	para[2] = FG_input[2]; // T_fluegas_in (K)
	para[3] = FG_input[3]; // Q_fluegas_in (kW)
	para[4] = vt[19] / 1.0e3; // Q_used (kW)
	para[5] = para[3] - para[4]; // Q_fluegas_out (kW)
	para[6] = para[5] / para[0]; // h_fluegas_out (kJ/kg)
	para[7] = Cal_T_FG(para[6]); // T_fluegas_out (K)

	st.insert(st.end(), LiBr_FG_name.begin(), LiBr_FG_name.end());
	vt.insert(vt.end(), para.begin(), para.end());

	return{ st,vt };
}