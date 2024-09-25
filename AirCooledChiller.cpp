#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

vector<string> ACC_name = {
	"Qeva,0", "COP0", "Qeva", "Teva,out",
	"Tcon,in", "¦¤Teva", "¦¤Tcon", "¦Ç", "cp,eva",
	"cp,con", "¦Ñeva", "¦Ñcon", "¦¤Teva,ref", "¦¤Tcon,ref",
	"Qcon,0", "Tcon,out", "Teva,in", "COPmax",
	"COP", "P", "Qcon", "Meva", "Mcon", "Veva", "Vcon"
	};

vector<vector<double>> ACC_Para_Set(const vector<double>& input) {
	
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
	para[9] = 1006.9; // cp,con (J/[kg¡¤K], under 41¡æ)
	para[10] = 999.74; // ¦Ñeva (kg/m3, under 9.5¡æ)
	para[11] = 1.1092; // ¦Ñcon (kg/m3, under 41¡æ)
	para[12] = 2.0; // ¦¤Teva,ref (K) modeling ?
	para[13] = 5.0; // ¦¤Tcon,ref (K) modeling ?
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

vector<double> ACC_AE_Fun(const vector<double>&x, const vector<double>& K) {
	vector<double> f(x.size()); // COPmax (-)
	f[0] = 1.0 / (((K[1] - K[2] * (1.0 / K[3] / x[0] + 1.0)) / K[0]) - 1.0) - x[0];

	return f;
}

vector<double> ACC_Var_Derive(const vector<double>& para, const vector<double>& sol) {
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

tuple<vector<string>, vector<double>> ACC_Sim(const vector<double>& input, bool start, vector<double> x0_pre) {
	vector<double> R;
	vector<vector<double>> Re1 = ACC_Para_Set(input);
	vector<vector<double>> Re2;
	vector<double> x0 = { 6.0 };

	if (start) {	
		Re2 = AE_Solve_AHS(ACC_AE_Fun, 20, x0, Re1[1], false);
	}
	else {
		Re2 = AE_Solve_Basic(ACC_AE_Fun, x0_pre, Re1[1]);
	}
	
	vector<double> Re3 = ACC_Var_Derive(Re1[0], Re2[0]);

	R.insert(R.end(), Re1[0].begin(), Re1[0].end());
	R.insert(R.end(), Re2[0].begin(), Re2[0].end());
	R.insert(R.end(), Re3.begin(), Re3.end());

	return { ACC_name,R };
}