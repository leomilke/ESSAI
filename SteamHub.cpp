#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_Load.h"

tuple<vector<string>, vector<double>> SteamHub_Sim(const vector<double>& T_in, const vector<double>& H_in, const vector<double>& M_in, const vector<double>& P_in) {
	int n = T_in.size();
	int len = 6 * n + 6;
	vector<double> Re;
	vector<double> X_in(n), S_in(n);
	vector<double> port_Mass_Flow_Rate = Scale(M_in, -1.0);
	Re.insert(Re.end(), T_in.begin(), T_in.end());
	Re.insert(Re.end(), H_in.begin(), H_in.end());
	Re.insert(Re.end(), port_Mass_Flow_Rate.begin(), port_Mass_Flow_Rate.end());
	Re.insert(Re.end(), P_in.begin(), P_in.end());
	
	for (double i = 0; i < n; ++i) { 
		if (T_in[i] == 0.0) {
			X_in[i] = 0.0;
			S_in[i] = 0.0;
		}
		else {
			X_in[i] = x_Ph(P_in[i], H_in[i]);
			S_in[i] = s_Px(P_in[i], X_in[i]);
		}
	}
	Re.insert(Re.end(), X_in.begin(), X_in.end());
	Re.insert(Re.end(), S_in.begin(), S_in.end());

	vector<double> Product = Multiply(H_in, M_in);
	vector<double> Product2 = Multiply(S_in, M_in);
	
	double M_SUM = accumulate(M_in.begin(), M_in.end(), 0.0);
	double h_OUT = accumulate(Product.begin(), Product.end(), 0.0) / M_SUM;
	double s_OUT = accumulate(Product2.begin(), Product2.end(), 0.0) / M_SUM * 1.03;

	Re.push_back(s_OUT);
	Re.push_back(-1.0);
	Re.push_back(h_OUT);
	Re.push_back(M_SUM);
	Re.push_back(T_hs(h_OUT, s_OUT));
	Re.push_back(Cal_P_Steam_T(*(Re.end() - 1.0)));
	Re[len - 5] = x_Ph(*(Re.end() - 1.0), h_OUT);

	vector<string> name(len);
	for (int i = 0; i < n; ++i) {
		name[i] = string("T_port_in_") + to_string(i + 1);
	}
	for (int i = n; i < 2 * n; ++i) {
		name[i] = string("H_port_in_") + to_string(i - n + 1);
	}
	for (int i = 2 * n; i < 3 * n; ++i) {
		name[i] = string("M_port_in_") + to_string(i - 2 * n + 1);
	}
	for (int i = 3 * n; i < 4 * n; ++i) {
		name[i] = string("P_port_in_") + to_string(i - 3 * n + 1);
	}
	for (int i = 4 * n; i < 5 * n; ++i) {
		name[i] = string("X_port_in_") + to_string(i - 4 * n + 1);
	}
	for (int i = 5 * n; i < 6 * n; ++i) {
		name[i] = string("S_port_in_") + to_string(i - 5 * n + 1);
	}
	name[len - 6] = "S_port_out";
	name[len - 5] = "X_port_out";
	name[len - 4] = "H_port_out";
	name[len - 3] = "M_port_out";
	name[len - 2] = "T_port_out";
	name[len - 1] = "P_port_out";

	return { name,Re };
}