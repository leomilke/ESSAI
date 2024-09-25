#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

extern deque<vector<double>> K_CoolingTower;

vector<string> CTowerClosed_name = {
	"M_coolingwater0_CT", "Q_air0_CT", "T_water_inlet0_CT", "T_water_out0_CT",
	"T_air_inlet0_CT", "T_air_out0_CT", "M_coolingwater_CT", "Q_air_CT",
	"T_water_inlet_CT", "T_air_inlet_CT", "M_spwater_loss0",
	"T_water_out_CT", "T_air_out_CT", "M_spwater_loss"
};

vector<double> CTowerClosed_Para_Set(const vector<double>& input) {
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
	para[10] = input[8]; // M_spwater_loss0 (kg/s)

	para[5] = para[4] + (4179.4 * para[0] * (para[2] - para[3]) - para[10] * 2257.0 * 1.0e3) / (para[1] / 3600.0 * 1.1495 * 1006.5); // T_air_out0_CT (K), water under 33¡æ, air under 30¡æ

	return para;
}

tuple<vector<string>, vector<double>> CTowerClosed_Simple_Sim(const vector<double>& input, bool start, bool iter_start) {
	vector<double> R, para = CTowerClosed_Para_Set(input), R1(3), T_w_out_st, e;
	double K = (para[2] - para[3]) / (para[5] - para[4]);

	for (double T_w_out = para[8]; T_w_out > para[9]; T_w_out = T_w_out - 0.1) {
		R1[2] = para[10] * (para[8] - T_w_out) / (para[2] - para[3]) * para[6] / para[0]; // M_spwater_loss (kg/s)
		R1[1] = (para[8] - T_w_out) / K + para[9]; // T_air_out_CT (K)
		e.push_back(abs(para[7] - ((4179.4 * para[6] * (para[8] - T_w_out) - R1[2] * 2257.0 * 1.0e3) / (R1[1] - para[9]) * 3600.0 / 1.1495 / 1006.5)));
		T_w_out_st.push_back(T_w_out);
	}

	R1[0] = T_w_out_st[distance(e.begin(), min_element(e.begin(), e.end()))];

	R.insert(R.end(), para.begin(), para.end());
	R.insert(R.end(), R1.begin(), R1.end());

	return { CTowerClosed_name,R };
}

vector<double> CTowerClosed_ODE_Fun(const vector<double>& x, double h, const vector<double>& K) {
	vector<double> dx(x.size());
	dx[0] = K[0] * (K[1] - x[0]); // Tw (K)
	
	return dx;
}

vector<double> CTowerClosed_AE_Fun(const vector<double>& K0, const vector<double>& para) {
	vector<double> f(1), K = K0;
	vector<vector<double>> R;
	vector<double> x0 = { para[2] };
	double range_CT[2] = { 0.0, 1.0 };
	K[1] = (para[4] + para[5]) / 2.0;

	R = ODE_Solve_Euler_Forward(CTowerClosed_ODE_Fun, range_CT, 0.01, x0, K);

	f[0] = *(R[1].end() - 1) - para[3];

	return f;
}

tuple<vector<string>, vector<double>> CTowerClosed_Sim(const vector<double>& input, bool start, bool iter_start) {
	vector<double> R, para = CTowerClosed_Para_Set(input), R1(3), T_air_out_st, e;
	double range_CT[2] = { 0.0, 1.0 };
	
	if (start && iter_start) {
		K_CoolingTower.push_front(AE_Solve_PSO(CTowerClosed_AE_Fun, 5.0, 2, 1e-2, para, 0.2)[0]);
		cout << "Cooling Tower Init completed !" << endl;
	}

	for (double T_air_out = para[9] + 0.1; T_air_out < para[8]; T_air_out = T_air_out + 0.1) {
		R1[0] = *(ODE_Solve_Euler_Forward(CTowerClosed_ODE_Fun, range_CT, 0.01, { para[8] }, { K_CoolingTower.front()[0], (T_air_out + para[9]) / 2.0 })[1].end() - 1); // T_water_out_CT (K)
		R1[2] = para[10] * (para[8] - R1[0]) / (para[2] - para[3]) * para[6] / para[0]; // M_spwater_loss (kg/s)
		e.push_back(abs(para[7] - ((4179.4 * para[6] * (para[8] - R1[0]) - R1[2] * 2257.0 * 1.0e3) / (T_air_out - para[9]) * 3600.0 / 1.1495 / 1006.5)));
		T_air_out_st.push_back(T_air_out);
	}

	R1[1] = T_air_out_st[distance(e.begin(), min_element(e.begin(), e.end()))];

	R.insert(R.end(), para.begin(), para.end());
	R.insert(R.end(), R1.begin(), R1.end());

	return { CTowerClosed_name,R };
}