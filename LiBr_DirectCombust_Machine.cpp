#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_Load.h"

vector<string> LiBr_DC_name = {
	"Q_gas", "V_gas_in", "¦Á", "Q_input", "Q_used",
	"Q_fluegas_out", "M_fluegas_out", "H_fluegas_out", "T_fluegas"
};

tuple<vector<string>, vector<double>> LiBr_DC_Sim(const vector<double>& input_mix, bool start, vector<double>& x0_pre, string mode) {
	tuple<vector<string>, vector<double>> temp;
	vector<double> para(9);
	vector<double> input;
	input.assign(input_mix.begin() + 3, input_mix.end() - 1);
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

	para[0] = input_mix[0]; // Q_gas (kJ/Nm3)
	para[1] = input_mix[1]; // V_gas_in (Nm3/h)
	para[2] = input_mix[2]; // ¦Á (-)
	para[3] = para[0] * para[1] / 3600.0; // Q_input (kW)
	para[4] = vt[19] / 1.0e3; // Q_used (kW)
	para[5] = para[3] - para[4]; // Q_fluegas_out (kW)
	para[6] = (para[1] * 0.7143 + para[2] * 9.5 * para[1] * 1.293) / 3600.0; // M_fluegas_out (kg/s)
	para[7] = para[5] / para[6]; // h_fluegas_out (kJ/kg)
	para[8] = Cal_T_FG(para[7]); // T_fluegas (K)

	st.insert(st.end(), LiBr_DC_name.begin(), LiBr_DC_name.end());
	vt.insert(vt.end(), para.begin(), para.end());

	return{ st,vt };
}