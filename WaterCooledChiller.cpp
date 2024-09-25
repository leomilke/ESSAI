#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

extern deque<vector<double>> K_CoolingTower;

vector<string> WCC_name = {
	"Qeva,0", "COP0", "Qeva", "Teva,out",
	"Tcon,in", "¦¤Teva", "¦¤Tcon", "¦Ç", "cp,eva",
	"cp,con", "¦Ñeva", "¦Ñcon", "¦¤Teva,ref", "¦¤Tcon,ref",
	"Qcon,0", "Tcon,out", "Teva,in", "COPmax",
	"COP", "P", "Qcon", "Meva", "Mcon", "Veva", "Vcon"
};

vector<vector<double>> WCC_Para_Set(const vector<double>& input) {
	
	vector<vector<double>> Re(2);
	vector<double> para(17);
	vector<double> K(4);
	
	para[0] = input[0]; // Qeva,0 (W)
	para[1] = input[1]; // COP0 (-)
	para[2] = input[2]; // Qeva (W)
	para[3] = input[3]; // Teva,out (K)
	para[4] = input[4]; // Tcon,in (K)
	para[5] = input[5]; // ¦¤Teva (K)
	para[6] = input[6]; // ¦¤Tcon (K)
	para[8] = 4196.0; // cp,eva (J/[kg¡¤K], under 9.5¡æ)
	para[9] = 4179.4; // cp,con (J/[kg¡¤K], under 33¡æ)
	para[10] = 999.74; // ¦Ñeva (kg/m3, under 9.5¡æ)
	para[11] = 994.70; // ¦Ñcon (kg/m3, under 33¡æ)
	para[12] = 2.0; // ¦¤Teva,ref (K) modeling ?
	para[13] = 2.0; // ¦¤Tcon,ref (K) modeling ?
	para[14] = - para[0] * (1.0 + para[1]) / para[1]; // Qcon,0 (W)
	para[15] = para[4] + para[6]; // Tcon,out (K)
	para[16] = para[3] - para[5]; // Teva,in (K)
	
	if (input.size() == 8) {
		para[7] = input[7]; // ¦Ç (-)
	}
	else
	{
		para[7] = para[1] / ((para[3] - para[12]) / (para[15] + para[13] - (para[3] - para[12])));
	}

	K[0] = para[3] - para[12] * para[2] / para[0];
	K[1] = para[15];
	K[2] = para[2] * para[13] / para[14];
	K[3] = para[7];

	Re[0] = para;
	Re[1] = K;

	return Re;
}

vector<double> WCC_AE_Fun(const vector<double>&x, const vector<double>& K) {
	vector<double> f(x.size()); // COPmax (-)
	f[0] = 1.0 / (((K[1] - K[2] * (1.0 / K[3] / x[0] + 1.0)) / K[0]) - 1.0) - x[0];

	return f;
}

vector<double> WCC_Var_Derive(const vector<double>& para, const vector<double>& sol) {
	vector<double> Re(7);
	Re[0] = para[7] * sol[0]; // COP (-)
	Re[1] = - para[2] / Re[0]; // P (W)
	Re[2] = Re[1] - para[2]; // Qcon (W)
	Re[3] = para[2] / para[8] / para[5]; // Meva (kg/s)
	Re[4] = Re[2] / para[9] / para[6]; // Mcon (kg/s)
	Re[5] = Re[3] / para[10]; // Veva (m3/s)
	Re[6] = Re[4] / para[11]; // Vcon (m3/s)

	return Re;
}

tuple<vector<string>, vector<double>> WCC_Sim(const vector<double>& input, bool start, vector<double> x0_pre) {
	if (input.size() == 7 || input.size() == 8) {
		vector<double> R;
		vector<vector<double>> Re1 = WCC_Para_Set(input);
		vector<vector<double>> Re2;
		vector<double> x0 = { 6.0 };

		if (start) {
			Re2 = AE_Solve_AHS(WCC_AE_Fun, 20, x0, Re1[1], false);
		}
		else {
			Re2 = AE_Solve_Basic(WCC_AE_Fun, x0_pre, Re1[1]);
		}

		vector<double> Re3 = WCC_Var_Derive(Re1[0], Re2[0]);

		R.insert(R.end(), Re1[0].begin(), Re1[0].end());
		R.insert(R.end(), Re2[0].begin(), Re2[0].end());
		R.insert(R.end(), Re3.begin(), Re3.end());

		return { WCC_name,R };
	}
	else if (input.size() == 17 || input.size() == 18) {
		vector<double> R, Re3, e, dt_st, x0 = { 6.0 };
		vector<double> input_WCC, input_CT_Open;
		vector<vector<double>> Re1, Re2;
		tuple<vector<string>, vector<double>> Re4;
		vector<string> WCC_com_name = WCC_name;
		double use_simple = *(input.end() - 1);
		
		if (input.size() == 17) {
			input_WCC.assign(input.begin(), input.begin() + 7);
			input_CT_Open.assign(input.begin() + 7, input.end() - 1);
		}
		else if (input.size() == 18) {
			input_WCC.assign(input.begin(), input.begin() + 8);
			input_CT_Open.assign(input.begin() + 8, input.end() - 1);
		}
		bool inner = true;
		for (double dt = 3.0; dt <= 12.0; dt = dt + 0.1) {
			input_WCC[6] = dt;
			Re1 = WCC_Para_Set(input_WCC);
			if (start) {
				Re2 = AE_Solve_AHS(WCC_AE_Fun, 20, x0, Re1[1], false);
			}
			else {
				Re2 = AE_Solve_Basic(WCC_AE_Fun, x0_pre, Re1[1]);
			}
			Re3 = WCC_Var_Derive(Re1[0], Re2[0]);
			
			input_CT_Open[5] = Re3[4];
			input_CT_Open[6] = Re1[0][15];
			
			if (use_simple > 0) {
				if (use_simple == 1) {
					Re4 = CTowerOpen_Simple_Sim(input_CT_Open, start, inner);
				}
				else if (use_simple == 2) {
					Re4 = CTowerClosed_Simple_Sim(input_CT_Open, start, inner);
				}
			}
			else if (use_simple == 0) {// very slow
				Re4 = CTowerOpen_Sim(input_CT_Open, start, inner);
			}
			else if (use_simple == -1) {// very slow
				Re4 = CTowerClosed_Sim(input_CT_Open, start, inner);
			}
			if (inner) { inner = false; }

			dt_st.push_back(dt);
			e.push_back(abs(get<1>(Re4)[11] - Re1[0][4]));
		}
		
		input_WCC[6] = dt_st[distance(e.begin(), min_element(e.begin(), e.end()))];
		
		Re1 = WCC_Para_Set(input_WCC);
		if (start) {
			Re2 = AE_Solve_AHS(WCC_AE_Fun, 20, x0, Re1[1], false);
		}
		else {
			Re2 = AE_Solve_Basic(WCC_AE_Fun, x0_pre, Re1[1]);
		}
		Re3 = WCC_Var_Derive(Re1[0], Re2[0]);
		
		input_CT_Open[5] = Re3[4];
		input_CT_Open[6] = Re1[0][15];

		if (use_simple > 0) {
			if (use_simple == 1) {
				Re4 = CTowerOpen_Simple_Sim(input_CT_Open, start, false);
			}
			else if (use_simple == 2) {
				Re4 = CTowerClosed_Simple_Sim(input_CT_Open, start, false);
			}
		}
		else if (use_simple == 0) {
			Re4 = CTowerOpen_Sim(input_CT_Open, start, false);
		}
		else if (use_simple == -1) {
			Re4 = CTowerClosed_Sim(input_CT_Open, start, false);
		}

		R.insert(R.end(), Re1[0].begin(), Re1[0].end());
		R.insert(R.end(), Re2[0].begin(), Re2[0].end());
		R.insert(R.end(), Re3.begin(), Re3.end());
		R.insert(R.end(), get<1>(Re4).begin(), get<1>(Re4).end());
		WCC_com_name.insert(WCC_com_name.end(), get<0>(Re4).begin(), get<0>(Re4).end());
		
		if (use_simple <= 0) {
			K_CoolingTower.push_back(K_CoolingTower.front());
			K_CoolingTower.pop_front();
		}
		
		return { WCC_com_name,R };
	}
}