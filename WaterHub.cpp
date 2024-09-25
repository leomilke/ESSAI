#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

tuple<vector<string>, vector<double>> WHub_Sim(const vector<double>& T_in, const vector<double>& M_in) {
	int n = T_in.size();
	int len = 2 * n + 2;
	vector<double> Re;
	vector<double> port_Mass_Flow_Rate = Scale(M_in, -1.0);
	Re.insert(Re.end(), T_in.begin(), T_in.end());
	Re.insert(Re.end(), port_Mass_Flow_Rate.begin(), port_Mass_Flow_Rate.end());

	vector<double> Product = Multiply(T_in, M_in);

	double M_SUM = accumulate(M_in.begin(), M_in.end(), 0.0);
	Re.push_back(accumulate(Product.begin(), Product.end(), 0.0) / M_SUM);
	Re.push_back(M_SUM);

	vector<string> name(len);
	for (int i = 0; i < n; ++i) { 
		name[i] = string("T_port_in_") + to_string(i + 1);
	}
	for (int i = n; i < 2*n; ++i) {
		name[i] = string("M_port_in_") + to_string(i - n + 1);
	}
	name[len - 2] = "T_port_out";
	name[len - 1] = "M_port_out";

	return { name,Re };
}