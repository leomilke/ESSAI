#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

tuple<vector<string>, vector<double>> GSHP_Sim(const vector<double>& input_constant, string mode, bool start, vector<double> x0_pre) {
	vector<double> input_0 = input_constant;
	vector<double> input_GHE;
	tuple<vector<string>, vector<double>> R, R_GHE;
	vector<double> input;
	vector<double> dt_st, error;
	double e, e_min_temp = 100.0, Tube_num = 0.0, dh_num = 0.0;
	if (input_0.size() == 14) {
		Tube_num = input_0[12];
		dh_num = input_0[13];
	}
	else if (input_0.size() == 15) {
		Tube_num = input_0[13];
		dh_num = input_0[14];
	}
	input_GHE.assign(input_0.end() - 7, input_0.end() - 2);
	input.assign(input_0.begin(), input_0.end() - 7);
	
	while (e_min_temp > 1.0) {
		for (double dt = 4.0; dt <= 16.0; dt = dt + 0.5) {
			if (mode == "Heating") {
				input[6] = -dt;
				R = WSHP_Sim(input, start, x0_pre);
			}
			else if (mode == "Cooling") {
				input[6] = dt;
				R = WCC_Sim(input, start, x0_pre);
			}
			R_GHE = GHE_Sim(input_GHE, get<1>(R)[22] / Tube_num, get<1>(R)[15], dh_num, mode);
			e = abs(get<1>(R_GHE)[17] - input[4]);
			error.push_back(e);
			dt_st.push_back(dt);
		}
		e_min_temp = *min_element(error.begin(), error.end());
		if (e_min_temp > 1.0) { 
			Tube_num = round(Tube_num / 2.0); 
			error.clear(); 
			dt_st.clear(); 
		}
	}
	
	input_0[input_0.size() - 2] = Tube_num;
	input_GHE.assign(input_0.end() - 7, input_0.end() - 2);

	if (mode == "Heating") {
		input[6] = -dt_st[distance(error.begin(), min_element(error.begin(), error.end()))];
		R = WSHP_Sim(input, start, x0_pre);
	}
	else if (mode == "Cooling") {
		input[6] = dt_st[distance(error.begin(), min_element(error.begin(), error.end()))];
		R = WCC_Sim(input, start, x0_pre);
	}
	R_GHE = GHE_Sim(input_GHE, get<1>(R)[22] / Tube_num, get<1>(R)[15], dh_num, mode);

	vector<string> OutS;
	vector<double> OutR;

	OutS.insert(OutS.end(), get<0>(R).begin(), get<0>(R).end());
	OutS.insert(OutS.end(), get<0>(R_GHE).begin(), get<0>(R_GHE).end());

	OutR.insert(OutR.end(), get<1>(R).begin(), get<1>(R).end());
	OutR.insert(OutR.end(), get<1>(R_GHE).begin(), get<1>(R_GHE).end());

	OutS.push_back("n_Tube");
	OutR.push_back(Tube_num);

	return { OutS,OutR };

}