#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

vector<string> WSHP_name = {
	"Qcon,0", "COP0", "Qcon", "Tcon,out",
	"Teva,in", "¦¤Tcon", "¦¤Teva", "¦Ç", "cp,con",
	"cp,eva", "¦Ñcon", "¦Ñeva", "¦¤Tcon,ref", "¦¤Teva,ref",
	"Qeva,0", "Teva,out", "Tcon,in", "COPmax",
	"COP", "P", "Qeva", "Mcon", "Meva", "Vcon", "Veva"
	};

vector<vector<double>> WSHP_Para_Set(const vector<double>& input) {
	
	vector<vector<double>> Re(2);
	vector<double> para(17);
	vector<double> K(4);
	
	para[0] = input[0]; // Qcon,0 (W)
	para[1] = input[1]; // COP0 (-)
	para[2] = input[2]; // Qcon (W)
	para[3] = input[3]; // Tcon,out (K)
	para[4] = input[4]; // Teva,in (K)
	para[5] = input[5]; // ¦¤Tcon (K)
	para[6] = input[6]; // ¦¤Teva (K)
	para[8] = 4179.4; // cp,con (J/[kg¡¤K], under 40¡æ)
	para[9] = 4195.2; // cp,eva (J/[kg¡¤K], under 10¡æ)
	para[10] = 992.22; // ¦Ñcon (kg/m3, under 40¡æ)
	para[11] = 999.70; // ¦Ñeva (kg/m3, under 10¡æ)
	para[12] = 2.0; // ¦¤Tcon,ref (K) modeling ?
	para[13] = 2.0; // ¦¤Teva,ref (K) modeling ?
	para[14] = para[0] * (1.0 - para[1]) / para[1]; // Qeva,0 (W)
	para[15] = para[4] + para[6]; // Teva,out (K)
	para[16] = para[3] - para[5]; // Tcon,in (K)
	
	if (input.size() == 8) {
		para[7] = input[7]; // ¦Ç (-)
	}
	else
	{
		para[7] = para[1] / ((para[3] + para[12]) / (para[3] + para[12] - (para[15] - para[13])));
	}

	K[0] = para[3] + para[12] * para[2] / para[0];
	K[1] = para[15];
	K[2] = para[2] * para[13] / -para[14];
	K[3] = para[7];

	Re[0] = para;
	Re[1] = K;

	return Re;
}

vector<double> WSHP_AE_Fun(const vector<double>&x, const vector<double>& K) {
	vector<double> f(x.size()); // COPmax (-)
	f[0] = 1.0 / (1.0 - ((K[1] + K[2] * (1.0 / K[3] / x[0] - 1.0)) / K[0])) - x[0];

	return f;
}

vector<double> WSHP_Var_Derive(const vector<double>& para, const vector<double>& sol) {
	vector<double> Re(7);
	Re[0] = para[7] * sol[0]; // COP (-)
	Re[1] = para[2] / Re[0]; // P (W)
	Re[2] = Re[1] - para[2]; // Qeva (W)
	Re[3] = para[2] / para[8] / para[5]; // Mcon (kg/s)
	Re[4] = Re[2] / para[9] / para[6]; // Meva (kg/s)
	Re[5] = Re[3] / para[10]; // Vcon (m3/s)
	Re[6] = Re[4] / para[11]; // Veva (m3/s)

	return Re;
}

tuple<vector<string>, vector<double>> WSHP_Sim(const vector<double>& input, bool start, vector<double> x0_pre) {
	vector<double> R;
	vector<vector<double>> Re1 = WSHP_Para_Set(input);
	vector<vector<double>> Re2;
	vector<double> x0 = { 6.0 };

	if (start) {
		Re2 = AE_Solve_AHS(WSHP_AE_Fun, 20, x0, Re1[1], false);
	}
	else {
		Re2 = AE_Solve_Basic(WSHP_AE_Fun, x0_pre, Re1[1]);
	}

	vector<double> Re3 = WSHP_Var_Derive(Re1[0], Re2[0]);
	
	R.insert(R.end(), Re1[0].begin(), Re1[0].end());
	R.insert(R.end(), Re2[0].begin(), Re2[0].end());
	R.insert(R.end(), Re3.begin(), Re3.end());

	return { WSHP_name,R };
}