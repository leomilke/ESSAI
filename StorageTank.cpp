#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

vector<string> STank_name = {
	"M_tank", "T_inlet_top/bottom", "H_tank", "A_tank", 
	"dh", "cp", "¦Ñ", "¦Ë", "u_in_tank", "T_outlet_top/bottom"
};

tuple<vector<string>, vector<double>, vector<double>> STank_Sim(const vector<double>& input, const vector<double>& T0, double dt, string mode) {
	vector<double> Re(10);
	Re[0] = input[0]; // M_tank (kg/s, >0 up  <0 down) 
	Re[1] = input[1]; // T_inlet_top/bottom (K)
	Re[2] = input[2]; // H_tank (m)
	Re[3] = input[3]; // A_tank (m2)
	double num = T0.size();
	Re[4] = Re[2] / (num - 1.0); // dh (m)
	if (mode == "Heating") {
		Re[5] = 4179.4; // cp (J/kg¡¤K, under 40¡æ)
		Re[6] = 992.22; // ¦Ñ (kg/m3, under 40¡æ)
		Re[7] = 0.62848; // ¦Ë (W/m¡¤K, under 40¡æ)
	}
	else if (mode == "Cooling") {
		Re[5] = 4196.0; // cp (J/kg¡¤K, under 9.5¡æ)
		Re[6] = 999.74; // ¦Ñ (kg/m3, under 9.5¡æ)
		Re[7] = 0.57772; // ¦Ë (W/m¡¤K, under 9.5¡æ)
	}
	Re[8] = Re[0] / Re[3] / Re[6]; // u_in_tank (m/s)
	
	double K1 = Re[8] * dt / Re[4];
	double K2 = Re[7] / Re[6] / Re[5] * dt / Re[4] / Re[4];

	vector<double> Tnew(num);
	vector<double> T0_use = T0;
	if (Re[0] > 0.0) {
		T0_use.push_back(T0[num - 1]);
		Tnew[0] = Re[1];
		for (int i = 1; i < num; ++i) {
			Tnew[i] = T0_use[i] + K1 * (T0_use[i - 1] - T0_use[i]) + K2 * (T0_use[i - 1] + T0_use[i + 1] - 2.0 * T0_use[i]);
		}
		Re[9] = Tnew[num - 1]; // T_outlet (K)
	}
	else if(Re[0] < 0.0) {
		Tnew[num - 1] = Re[1];
		for (int i = num - 2; i > 0; --i) {
			Tnew[i] = T0_use[i] + K1 * (T0_use[i] - T0_use[i + 1]) + K2 * (T0_use[i - 1] + T0_use[i + 1] - 2.0 * T0_use[i]);
		}
		Tnew[0] = T0_use[0] + K1 * (T0_use[0] - T0_use[1]) + K2 * (T0_use[0] + T0_use[1] - 2.0 * T0_use[0]);
		Re[9] = Tnew[0]; // T_outlet (K)
	}
	else if (Re[0] == 0.0) {
		T0_use.push_back(T0[num - 1]);
		Tnew[0] = T0_use[0] + K2 * (T0_use[0] + T0_use[1] - 2.0 * T0_use[0]);
		for (int i = 1; i < num; ++i) {
			Tnew[i] = T0_use[i] + K2 * (T0_use[i - 1] + T0_use[i + 1] - 2.0 * T0_use[i]);
		}
		Re[9] = Tnew[num - 1]; // T_outlet (K)
	}

	return { STank_name,Re,Tnew };
}