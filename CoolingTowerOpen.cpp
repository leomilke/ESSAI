#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

extern deque<vector<double>> K_CoolingTower;

vector<string> CTowerOpen_name = {
	"M_coolingwater0_CT", "Q_air0_CT", "T_water_inlet0_CT", "T_water_out0_CT",
	"T_air_inlet0_CT", "T_air_out0_CT", "M_coolingwater_CT", "Q_air_CT",
	"T_water_inlet_CT", "T_air_inlet_CT", "R_water_drift0",
	"T_water_out_CT", "T_air_out_CT", "R_water_drift"
};

vector<double> CTowerOpen_Para_Set(const vector<double>& input) {
	vector<double> para(11);
	para[0] = input[0]; // M_coolingwater0_CT (kg/s)
	para[1] = input[1]; // Q_air0_CT (m3/h)
	para[2] = input[2]; // T_water_inlet0_CT (K)
	para[3] = input[3]; // T_water_out0_CT (K)
	para[4] = input[4]; // T_air_inlet0_CT (K)
	
	para[6] = input[5]; // M_coolingwater_CT (kg/s)
	para[7] = para[1] / para[0] * para[6]; // Q_air_CT (m3/h), not final value
	para[8] = input[6]; // T_water_inlet_CT (K)
	para[9] = input[7]; // T_air_inlet_CT (K)
	para[10] = input[8]; // R_water_drift0 (%)
	
	para[5] = para[4] + (4179.4 * para[0] * (para[2] - para[3]) - para[10] / 100.0 * para[0] * 2257.0 * 1.0e3) / (para[1] / 3600.0 * 1.1495 * 1006.5); // T_air_out0_CT (K), water under 33¡æ, air under 30¡æ
	
	return para;
}

vector<double> CTowerOpen_ODE_Fun(const vector<double>& x, double h, const vector<double>& K) {
	vector<double> dx(x.size());
	double Pair = x[2] * 100000.0 / 217.0 * x[1]; // Pa
	double Psat = 611.2 * exp((17.62 * (x[1] - 273.15)) / (243.12 + (x[1] - 273.15))); // Pa

	dx[0] = -K[0] * (x[0] - x[1]) - K[1] / 1.0e4 * (Psat - Pair); // Tw (K)
	dx[1] = -K[2] * (x[0] - x[1]); // Ta (K)
	dx[2] = -K[3] / 1.0e6 * (Psat - Pair); // m (kg/m3)
	
	return dx;
}

vector<double> CTowerOpen_AE_Fun(const vector<double>& K, const vector<double>& para) {
	vector<double> f(6);
	vector<vector<double>> R;
	vector<double> x0 = { para[2], para[5], para[10] / 100.0 * para[0] / (para[1] / 3600.0) };
	double range_CT[2] = { 0.0, 1.0 };

	R = ODE_Solve_Euler_Forward(CTowerOpen_ODE_Fun, range_CT, 0.01, x0, K);

	f[0] = *(R[1].end() - 1) - para[3];
	f[1] = *(R[2].end() - 1) - para[4];

	x0 = { para[2], para[5], 0.5 * para[10] / 100.0 * para[0] / (para[1] / 3600.0) };

	R = ODE_Solve_Euler_Forward(CTowerOpen_ODE_Fun, range_CT, 0.01, x0, K);

	f[2] = *(R[1].end() - 1) - para[3];
	f[3] = *(R[2].end() - 1) - para[4];

	x0 = { para[2], para[5], 0.2 * para[10] / 100.0 * para[0] / (para[1] / 3600.0) };

	R = ODE_Solve_Euler_Forward(CTowerOpen_ODE_Fun, range_CT, 0.01, x0, K);

	f[4] = *(R[1].end() - 1) - para[3];
	f[5] = *(R[2].end() - 1) - para[4];

	return f;
}

tuple<vector<string>, vector<double>> CTowerOpen_Sim(const vector<double>& input, bool start, bool iter_start) {
	vector<double> R, Re3(3), Re1 = CTowerOpen_Para_Set(input), x0;
	vector<vector<double>> Para_inside, Para_inside_op;
	double range_CT[2] = { 0.0, 1.0 };
	double emin = 1.0e10, et, T_air_out_op = 0.0, T_w_out = Re1[3], Q_air_CT_cal, Rd = Re1[10] * (Re1[8] - T_w_out) / (Re1[2] - Re1[3]);
	
	if (start && iter_start) {
		K_CoolingTower.push_front(AE_Solve_PSO(CTowerOpen_AE_Fun, 4.0, 4, 1.0e-1, Re1, 0.01)[0]);
		cout << "Cooling Tower Init completed !" << endl;
	}

	do {
		Q_air_CT_cal = Re1[7];
		x0 = { Re1[8], 0.0, Rd / 100.0 * Re1[6] / (Re1[7] / 3600.0) };
		for (double T_ao = Re1[9]; T_ao <= Re1[8]; T_ao = T_ao + 0.1) {
			x0[1] = T_ao;
			Para_inside = ODE_Solve_Euler_Forward(CTowerOpen_ODE_Fun, range_CT, 0.01, x0, K_CoolingTower.front());
			et = abs(*(Para_inside[2].end() - 1) - Re1[9]);
			if (et < emin) {
				emin = et;
				T_air_out_op = T_ao;
				T_w_out = *(Para_inside[1].end() - 1);
				Para_inside_op = Para_inside;
			}
		}
		Rd = Re1[10] * (Re1[8] - T_w_out) / (Re1[2] - Re1[3]);
		Re1[7] = (4179.4 * Re1[6] * (Re1[8] - T_w_out) - Rd / 100.0 * Re1[6] * 2257.0 * 1.0e3) / (T_air_out_op - Re1[9]) * 3600.0 / 1.1495 / 1006.5;
	} while (abs(Re1[7] - Q_air_CT_cal) > 1.0);
	
	Re3[0] = T_w_out; // T_water_out_CT (K)
	Re3[1] = T_air_out_op; // T_air_out_CT (K)
	Re3[2] = Rd; // R_water_drift (%)

	R.insert(R.end(), Re1.begin(), Re1.end());
	R.insert(R.end(), Re3.begin(), Re3.end());

	return { CTowerOpen_name,R };
}

tuple<vector<string>, vector<double>> CTowerOpen_Simple_Sim(const vector<double>& input, bool start, bool iter_start) {
	vector<double> R, para = CTowerOpen_Para_Set(input), R1(3), T_w_out_st, e;
	double K = (para[2] - para[3]) / (para[5] - para[4]);
	
	for (double T_w_out = para[8]; T_w_out > para[9]; T_w_out = T_w_out - 0.1) {
		R1[2] = para[10] * (para[8] - T_w_out) / (para[2] - para[3]); // R_water_drift (%)
		R1[1] = (para[8] - T_w_out) / K + para[9]; // T_air_out_CT (K)
		e.push_back(abs(para[7] - ((4179.4 * para[6] * (para[8] - T_w_out) - R1[2] / 100.0 * para[6] * 2257.0 * 1.0e3) / (R1[1] - para[9]) * 3600.0 / 1.1495 / 1006.5)));
		T_w_out_st.push_back(T_w_out);
	}
	
	R1[0] = T_w_out_st[distance(e.begin(), min_element(e.begin(), e.end()))];

	R.insert(R.end(), para.begin(), para.end());
	R.insert(R.end(), R1.begin(), R1.end());
	
	return { CTowerOpen_name,R };
}