#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"

vector<string> GHE_name = {
	"Rtube,out", "Rtube,in", "Htube", "¦Ëtube",
	"Mf_GHE", "T_GHE_in", "Tg", "¦Ñ_f_GHE", "cp_f_GHE",
	"¦Ë_f_GHE", "¦Ì_f_GHE", "v_f_GHE", "vel_f_GHE",
	"Re number", "Pr number", "Nu number", "h_f_in",
	"T_GHE_out"
};

tuple<vector<string>, vector<double>> GHE_Sim(const vector<double>& input, double M, double T_in, double dh_num, string mode) {
	vector<double> Re(18);
	double Ltube, Ain, PI=3.1415926, f, K, Tout=T_in;

	Re[0] = input[0]; // Rtube,out (m)
	Re[1] = input[1]; // Rtube,in (m)
	Re[2] = input[2]; // Htube (m)
	Re[3] = input[3]; // ¦Ëtube (W/m¡¤K)
	Re[4] = M; // Mf_GHE (kg/s)
	Re[5] = T_in; // T_GHE_in (K)
	Re[6] = input[4]; // Tg (K)

	if (mode == "Heating") {
		// under 10¡æ
		Re[7] = 999.70; // ¦Ñ_f_GHE (kg/m3)
		Re[8] = 4195.2; // cp_f_GHE (J/kg¡¤K)
		Re[9] = 0.57878; // ¦Ë_f_GHE (W/m¡¤K)
		Re[10] = 0.0013059; // ¦Ì_f_GHE (Pa¡¤s)
	}
	else if (mode == "Cooling") {
		// under 24¡æ
		Re[7] = 997.30; // ¦Ñ_f_GHE (kg/m3)
		Re[8] = 4181.8; // cp_f_GHE (J/kg¡¤K)
		Re[9] = 0.60487; // ¦Ë_f_GHE (W/m¡¤K)
		Re[10] = 0.00091068; // ¦Ì_f_GHE (Pa¡¤s)
	}
	Re[11] = Re[10] / Re[7]; // v_f_GHE (m^2/s)

	Ltube = 2 * Re[2];
	Ain = PI * Re[1] * Re[1];
	Re[12] = Re[4] / Ain / Re[7]; // vel_f_GHE (m/s)
	Re[13] = Re[12] * 2.0 * Re[1] / Re[11]; // Re number (-)
	Re[14] = Re[10] * Re[8] / Re[9]; // Pr number (-)
	
	if (Re[13] <= 2300) {
		Re[15] = 3.66; // Nu number (-) for laminar flow
	}
	else {
		f = pow((1.82 * log10(Re[13]) - 1.5), -2.0);
		Re[15] = f / 8.0 * (Re[13] - 1000.0) * Re[14] * (1.0 + pow(2.0 * Re[1] / Ltube, 2.0 / 3.0)) / (1.0 + 12.7 * sqrt(f / 8.0) * (pow(Re[14], 2.0 / 3.0) - 1.0)); // Gnielinski formula for turbulent flow
	}

	Re[16] = Re[15] / 2.0 / Re[1] * Re[9]; // h_f_in (W/m2¡¤K)

	double dh = Ltube / dh_num;
	K = 2 * PI * dh / M / Re[8] / (log(Re[0] / Re[1]) / Re[3] + 1.0 / Re[16] / Re[1]);
	for (double hh = 0.0; hh <= Ltube; hh = hh + dh) {
		Tout = (Tout + K * Re[6]) / (1 + K);
	}

	Re[17] = Tout; // T_GHE_out (K)

	return { GHE_name,Re };
}