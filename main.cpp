#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Sub_Models.h"
#include "Pre_Load.h"
#include <time.h>
using namespace std;

int main()
{
	double dt = 0.1, T = 4000; //86400
	
	vector<string> dev_HVAC = { "ASHP", "WSHP", "GSHP", "LiBr_DC", "CB_Boiler", "CFB_Boiler", "Gas_Boiler", "Elec_Boiler", "LiBr_WHW", "LiBr_FG", "Gas_Boiler", "LiBr_Solar" };
	vector<vector<double>> input_HVAC = {
		{520 * 1e3, 5.4, 500 * 1e3, 50.0 + 273.15, 5.0 + 273.15, 5.0, -10.0},
		{500 * 1e3, 4.5, 500 * 1e3, 45.0 + 273.15, 12.0 + 273.15, 5.0, -5.0},
		{580 * 1e3, 5.2, 560 * 1e3, 48.0 + 273.15, 15.0 + 273.15, 5.0, -5.0, 0.333, 0.03, 0.02, 100.0, 0.5, 20.0 + 273.15, 40.0, 1000.0},
		{35590.0, 50.0, 1.2, 520 * 1e3, 1.2, 500 * 1e3, 50.0 + 273.15, 12.0 + 273.15, 5.0, -6.0, 2},
		{29307.6, 0.015, 1.2, 20.0 + 273.15, 800.0 + 273.15, 300.0 + 273.15, 0.93, 0.95, 0.91, 0},
		{29307.6, 0.015, 1.2, 20.0 + 273.15, 800.0 + 273.15, 700.0 + 273.15, 0.93, 0.8, 0.85, 0},
		{35590.0, 10.0, 4, 20.0 + 273.15, 800.0 + 273.15, 300.0 + 273.15, 0.93, 0.9, 0.90, 1.2, 1.1, 0},
		{1 * 1e3, 0.9, 4.0, 20.0 + 273.15, 0.95},
		{100.0 + 273.15, 3.0, 500 * 1e3, 1.3, 500 * 1e3, 45.0 + 273.15, 12.0 + 273.15, 5.0, -5.0, 2},
		{530 * 1e3, 6.1, 500 * 1e3, 48.0 + 273.15, 5.0 + 273.15, 5.0, -12.0, 0.333},
		{35590.0, 10.0, 4, 22.0 + 273.15, 800.0 + 273.15, 300.0 + 273.15, 0.95, 0.85, 0.90, 1.2, 1.1, 0},
		{300.0, 1.8, -1, -1, -1, -1, 3.0, 10.0, 1000.0, 60.0 + 273.15, 70.0 + 273.15, 3.0, 520 * 1e3, 1.4, 500 * 1e3, 50.0 + 273.15, 5.0 + 273.15, 5.0, -10.0, 0 }
	};
	vector<double> input_ST = { 5.0, 50.0 + 273.15, 20.0, 5.0 };
	vector<double> T0(1000.0, 12.0 + 273.15);
	vector<int> charge = { 0, 1, 1, 0 ,1, 0, 0, 0, 0, 1, 0, 0 };
	vector<string> LB_S_TP = { };
	vector<string> LBFG_TP = { "Air" };
	vector<tuple<string, vector<int>>> LBFG_C{
		{ "Water", {5} }
	};
	vector<vector<double>> FGWasteHeatUse_W = {
		{ 4, 6 }, // WaterBoiler
		{ 2 }  // SteamBoiler
	};
	vector<vector<double>> input_WHHWB = {
		{0, 0.70, 0.03, 20.0 + 273.15, 0.0},
		{0, 0.65, 0.05, 15.0 + 273.15, 1.0},
		{0, 0.75, 0.04, 18.0 + 273.15, 0.0}
	};
	vector<vector<double>> FGWasteHeatUse_S = {
		{  }, // WaterBoiler
		{ 0, 1, 3 }  // SteamBoiler
	};
	vector<vector<double>> input_WHSB = {
		{0, 0.80, 0.015, 30.0 + 273.15, 0.65},
		{0, 0.75, 0.02, 25.0 + 273.15, 0.65},
		{0, 0.85, 0.01, 28.0 + 273.15, 0.6}
	};
	
	vector<string> dev_Steam = { "Gas_Steam_Boiler", "Gas_Steam_Boiler", "CFB_Steam_Boiler", "CB_Steam_Boiler" };
	vector<vector<double>> input_Steam = {
		{ 35590.0, 60.0, 0.8, 20.0 + 273.15, 800.0 + 273.15, 400.0 + 273.15, 0.93, 1.0, 0.93, 1.2, 1.2, 0.6, 0.6, 0.0 },
		{ 35590.0, 15.0, 0.2, 20.0 + 273.15, 800.0 + 273.15, 300.0 + 273.15, 0.93, 0.9, 0.90, 1.2, 1.1, 0.7, 0.65, 0.0 },
		{ 29307.6, 0.015, 0.12, 20 + 273.15, 800.0 + 273.15, 700.0 + 273.15, 0.93, 0.95, 0.91, 0.7, 0.7, 0.0 },
		{ 29307.6, 0.015, 0.12, 20 + 273.15, 800.0 + 273.15, 300.0 + 273.15, 0.93, 0.95, 0.91, 0.7, 0.7, 0.0 }
	};

	vector<vector<int>> LB_St = {};
	vector<string> LB_TP = {};

	clock_t s = clock();
	auto R = HVAC_Steam_Sim(T, dt, dev_HVAC, dev_Steam, input_HVAC, input_Steam, true, input_ST, T0, charge, false, 0, LB_TP, LB_St, LBFG_C, LBFG_TP, FGWasteHeatUse_W, input_WHHWB, FGWasteHeatUse_S, input_WHSB);
	clock_t e = clock();
	cout << e - s << endl;

	for (int i = 0; i < get<0>(get<1>(R)[10000]).size(); ++i) {
		cout << get<0>(get<1>(R)[10000])[i] << "    " << get<1>(get<1>(R)[10000])[i] << endl;
	}

	for (int i = 0; i < get<0>(get<2>(R)[1000]).size(); ++i) {
		cout << get<0>(get<2>(R)[1000])[i] << "    " << get<1>(get<2>(R)[1000])[i] << endl;
	}

	
	
	
	return 0;
}