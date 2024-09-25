#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_Load.h"

tuple<vector<string>, vector<double>> FGHub_Sim(const vector<double>& Q_in, const vector<double>& T_in) {
	int n = Q_in.size();
	int len = 4 * n + 4;
	vector<double> Re;
	vector<double> H_in(n), M_in(n);
	
	for (int i = 0; i < n; ++i) {
		if (Q_in[i] == 0.0) {
			H_in[i] = 0.0;
			M_in[i] = 0.0;
		}
		else {
			H_in[i] = Cal_h_FG(T_in[i]); // kJ/kg
			M_in[i] = -1.0 * Q_in[i] / H_in[i]; // kg/s 
		}
	}
	Re.insert(Re.end(), Q_in.begin(), Q_in.end());
	Re.insert(Re.end(), T_in.begin(), T_in.end());
	Re.insert(Re.end(), H_in.begin(), H_in.end());
	Re.insert(Re.end(), M_in.begin(), M_in.end());	
	
	double M_SUM = -accumulate(M_in.begin(), M_in.end(), 0.0);
	vector<double> Product = Multiply(H_in, M_in);
	double h_OUT = -accumulate(Product.begin(), Product.end(), 0.0) / M_SUM;
	double T_OUT = Cal_T_FG(h_OUT);
	double Q_OUT = h_OUT * M_SUM;
	
	Re.push_back(M_SUM);
	Re.push_back(h_OUT);
	Re.push_back(T_OUT);
	Re.push_back(Q_OUT);
	
	vector<string> name(len);
	for (int i = 0; i < n; ++i) {
		name[i] = string("Q_port_in_") + to_string(i + 1); // kW
	}
	for (int i = n; i < 2 * n; ++i) {
		name[i] = string("T_port_in_") + to_string(i - n + 1);
	}
	for (int i = 2 * n; i < 3 * n; ++i) {
		name[i] = string("H_port_in_") + to_string(i - 2 * n + 1);
	}
	for (int i = 3 * n; i < 4 * n; ++i) {
		name[i] = string("M_port_in_") + to_string(i - 3 * n + 1);
	}
	
	name[len - 4] = "M_port_out";
	name[len - 3] = "h_port_out";
	name[len - 2] = "T_port_out";
	name[len - 1] = "Q_port_out";

	return { name,Re };
}