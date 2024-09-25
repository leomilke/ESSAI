#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_Load.h"

extern deque<double> Load_pre;
extern deque<vector<double>> K_combust_boiler;

vector<string> GSB_name = {
	"Q_gas_fuel", "V_gas_fuel0", "M_steam", "T_water_in", "h_steam_out0",
	"T_boiler_bottom_0", "T_boiler_out_0", "eta_boiler_0", "Load",
	"V_gas_fuel", "eta_boiler", "h_steam_out", "¦Á0", "¦Á", "p_steam_out0", 
	"p_steam_out", "T_steam_out0", "T_steam_out", "T_boiler_bottom", "T_boiler_out"
};

vector<double> GSB_Para_Set(const vector<double>& input) {
	vector<double> para(18);
	para[0] = input[0]; // Q_gas_fuel (kJ/Nm3)
	para[1] = input[1]; // V_gas_fuel0 (Nm3/h)
	para[2] = input[2]; // M_steam (t/h)
	para[3] = input[3]; // T_water_in (K)
	para[5] = input[4]; // T_boiler_bottom_0 (K)
	para[6] = input[5]; // T_boiler_out_0 (K)
	para[7] = input[6]; // eta_boiler_0 (-)

	para[8] = input[7]; // Load (-)
	para[9] = para[8] * para[1]; // V_gas_fuel (Nm3/h)
	para[10] = input[8]; // eta_boiler (-)

	para[4] = 1.0e-3 * para[0] * para[1] * para[7] / para[2] + (para[3]-273.15) * 4.1841; // h_steam_out0 (kJ/kg), Cp_water under 20¡æ
	para[11] = 1.0e-3 * para[0] * para[9] * para[10] / para[2] + (para[3]-273.15) * 4.1841; // h_steam_out (kJ/kg)

	para[12] = input[9]; // ¦Á0 (-)
	para[13] = input[10]; // ¦Á (-)
	para[14] = input[11]; // p_steam_out0 (MPa)
	para[15] = input[12]; // p_steam_out (MPa)

	para[16] = Cal_T_Steam_P(para[14]); // T_steam_out0 (K)
	para[17] = Cal_T_Steam_P(para[15]); // T_steam_out (K)

	return para;
}

vector<double> GSB_ODE_Fun(const vector<double>& x, double h, const vector<double>& K) {
	vector<double> dx(x.size());
	double E = 0.01305e4, R = 8.314;
	dx[0] = K[0] * 650.0 * 5.0e2 * exp(-E / R / x[0]) * x[2] * x[3] - K[1] * (650.0 / 400.0) * (x[0] - x[1]); // T
	dx[1] = K[2] * (150.0 / 400.0) * (x[0] - x[1]); // Tst
	dx[2] = -2.0 * 4.0e3 * exp(-E / R / x[0]) * x[2] * x[3]; // fO2
	dx[3] = -4.0e3 * exp(-E / R / x[0]) * x[2] * x[3]; // fCH4

	return dx;
}

vector<double> GSB_AE_Fun(const vector<double>& K, const vector<double>& para) {
	vector<double> f(K.size());
	vector<vector<double>> R;
	vector<double> x0 = { para[5], para[3], 0.21 * para[12] / (para[12] + 0.105), 0.105 / (para[12] + 0.105) };
	double range_boiler[2] = { 0.0, 1.0 };

	R = ODE_Solve_Euler_Forward(GSB_ODE_Fun, range_boiler, 0.001, x0, K);

	f[0] = *(R[1].end() - 1) - para[6]; 
	f[1] = *(R[2].end() - 1) - para[16];
	if (any_of(R[3].begin(), R[3].end(), [](double x) {return x < 0; }) || any_of(R[4].begin(), R[4].end(), [](double x) {return x < 0; }))
	{
		f[2] = 1.0e9;
		f[3] = 1.0e9;
	}
	else {
		f[2] = (*(R[3].end() - 1) - 0.21 * (para[12] - 1.0) / (para[12] + 0.105)) * 1.0e3;
		f[3] = *(R[4].end() - 1) * 1.0e3;
	}
	
	return f;
}

tuple<vector<string>, vector<double>> GSB_Sim(const vector<double>& input, bool start, vector<double> T_b_b0) {
	vector<double> R;
	vector<double> Re1 = GSB_Para_Set(input);
	vector<double> x0 = { Re1[5], Re1[3], 0.21 * Re1[13] / (Re1[13] + 0.105), 0.105 / (Re1[13] + 0.105) };
	double range_boiler[2] = { 0.0, 1.0 };
	vector<vector<double>> Para_inside, Para_inside_op;
	double emin = 1.0e10, et, T_b_b_op = 0.0, T_b_o = 0.0;
	vector<double> Re3(2);

	if (start) {
		K_combust_boiler.push_front(AE_Solve_PSO(GSB_AE_Fun, 5.0, 4, 1.0e-3, Re1, 0.01)[0]);
		Load_pre.push_front(-1.0);
		cout << "Gas Steam Boiler Init completed !" << endl;
	}
	if (Re1[8] == Load_pre.front()) {
		T_b_b_op = T_b_b0[0];
		x0[0] = T_b_b_op;
		Para_inside_op = ODE_Solve_Euler_Forward(GSB_ODE_Fun, range_boiler, 0.001, x0, K_combust_boiler.front());
		T_b_o = *(Para_inside_op[1].end() - 1);
	}
	else {
		for (double T_b_b = 400.0 + 273.15; T_b_b <= Re1[5]; ++T_b_b) {
			x0[0] = T_b_b;
			Para_inside = ODE_Solve_Euler_Forward(GSB_ODE_Fun, range_boiler, 0.001, x0, K_combust_boiler.front());
			et = abs(*(Para_inside[2].end() - 1) - Re1[17]);
			if (et < emin) {
				emin = et;
				T_b_b_op = T_b_b;
				T_b_o = *(Para_inside[1].end() - 1);
				Para_inside_op = Para_inside;
			}
		}
	}
	K_combust_boiler.push_back(K_combust_boiler.front());
	K_combust_boiler.pop_front();
	Load_pre.push_back(Re1[8]);
	Load_pre.pop_front();
	Re3[0] = T_b_b_op; // T_boiler_bottom (K)
	Re3[1] = T_b_o; // T_boiler_out (K)

	R.insert(R.end(), Re1.begin(), Re1.end());
	R.insert(R.end(), Re3.begin(), Re3.end());

	return { GSB_name,R };
}

vector<string> GSB_simple_name = {
	"Q_fuel", "V_fuel0", "Load", "V_fuel", "M_steam", "T_water_in", 
	"h_steam_out", "eta_boiler", "p_steam_out", "T_steam_out"
};

tuple<vector<string>, vector<double>> GSB_Simple_Sim(const vector<double>& input) {
	vector<double> para(10);
	para[0] = input[0]; // Q_fuel (kJ/Nm3)
	para[1] = input[1]; // V_fuel0 (Nm3/h)
	para[2] = input[2]; // Load (-)
	para[3] = para[1] * para[2]; // V_fuel (Nm3/h)
	para[4] = input[3]; // M_steam (t/h)
	para[5] = input[4]; // T_water_in (K)
	para[7] = input[5]; // eta_boiler (-)

	para[6] = 1.0e-3 * para[0] * para[3] * para[7] / para[4] + (para[5]-273.15) * 4.1841; // h_steam_out (kJ/kg), Cp_water under 20¡æ	
	para[8] = input[6]; // p_steam_out (MPa)
	para[9] = Cal_T_Steam_P(para[8]); // T_steam_out (K)

	return { GSB_simple_name,para };
}