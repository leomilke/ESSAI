#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
// Better use small time step.
extern deque<vector<double>> ST_T_inside;

vector<string> LiBr_Solar_name = { "P_Aux_heat" };

tuple<vector<string>, vector<double>> LiBr_Solar_Sim(const vector<double>& input_mix, bool start, vector<double>& x0_pre, string mode, double dt) {
	tuple<vector<string>, vector<double>> temp_SWH, temp;
	tuple<vector<string>, vector<double>, vector<double>> temp_ST, temp_ST1;
	vector<double> input_SWH, input_ST, input_ST_use, input;
	input_SWH.assign(input_mix.begin(), input_mix.begin() + 4);
	input_ST.assign(input_mix.begin() + 4, input_mix.begin() + 10);
	input.assign(input_mix.begin() + 10, input_mix.end());
	double P_aux = 0.0;

	if (start) {
		vector<double> ST_T0(input_ST[4], input_ST[5]);
		ST_T_inside.push_front(ST_T0);
	}
	
	input_SWH[2] = ST_T_inside.front()[10];
	input_SWH[3] = input[1];
	temp_SWH = SWH_Sim(input_SWH, input[0]);
	
	input_ST[0] = -get<1>(temp_SWH)[4];
	input_ST[1] = get<1>(temp_SWH)[5];
	input_ST_use.assign(input_ST.begin(), input_ST.end() - 2);

	temp_ST = STank_Sim(input_ST_use, ST_T_inside.front(), dt, "Heating");
	
	input[0] = *(ST_T_inside.front().end() - 1 - 10);
	if (input[0] < 70.0 + 273.15) {
		P_aux = ((70.0 + 273.15) - input[0]) * 4.1850 * input[1]; // kW, cp 60¡æ
		input[0] = 70.0 + 273.15;
	}

	temp = LiBr_WHW_Sim(input, start, x0_pre, mode);
	
	input_ST_use[0] = *(get<1>(temp).end() - 5);
	input_ST_use[1] = *(get<1>(temp).end() - 1);

	temp_ST1 = STank_Sim(input_ST_use, get<2>(temp_ST), dt, "Heating");
	
	vector<string>st = get<0>(temp);
	vector<double>vt = get<1>(temp);
	
	st.insert(st.end(), get<0>(temp_SWH).begin(), get<0>(temp_SWH).end());
	vt.insert(vt.end(), get<1>(temp_SWH).begin(), get<1>(temp_SWH).end());
	st.insert(st.end(), get<0>(temp_ST).begin(), get<0>(temp_ST).end());
	vt.insert(vt.end(), get<1>(temp_ST).begin(), get<1>(temp_ST).end());
	st.insert(st.end(), get<0>(temp_ST1).begin(), get<0>(temp_ST1).end());
	vt.insert(vt.end(), get<1>(temp_ST1).begin(), get<1>(temp_ST1).end());
	st.insert(st.end(), LiBr_Solar_name.begin(), LiBr_Solar_name.end());
	vt.push_back(P_aux);

	ST_T_inside.push_back(get<2>(temp_ST1));
	ST_T_inside.pop_front();

	return{ st,vt };
}