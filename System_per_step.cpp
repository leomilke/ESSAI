#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_Load.h"

extern vector<vector<double>> FGInput_Steam, FGInput_Water;
extern vector<tuple<string, vector<int>>> LB_FG_Source;
extern deque<string> LB_FG_TYPE;
extern deque<int> LB_FG_INDEX;

tuple<vector<string>, vector<double>, vector<vector<double>>, vector<vector<double>>> SteamSysSim_single_step(const vector<string>& heater, const vector<vector<double>>& input, bool start, vector<vector<double>> x0_pre, map<int, int>& LB_St_Con, int LB_num, vector<vector<double>> FGWasteHeatUse_S, vector<vector<double>> input_WHSB) {

	vector<string> para_name;
	vector<double> value;
	int len = heater.size();
	vector<double> T_in(len), M_in(len), P_in(len), H_in(len);

	tuple<vector<string>, vector<double>> temp, temp_WHB;
	vector<string> st, vs_WHB, vs_WHB_temp;
	vector<double> vt, vt_FG, vt_WHB;
	vector<vector<double>> x0_pre_all;
	vector<vector<double>> Steam2HVAC(LB_num);
	vector<vector<vector<double>>> LB_Steam_part(LB_num);

	for (int i = 0; i < len; ++i) {
		vector<double> x0_pre_i;

		if (heater[i] == "Elec_Steam_Boiler") {
			temp = ESB_Sim(input[i]);
		}
		else if (heater[i] == "CB_Steam_Boiler") {
			if (*(input[i].end() - 1) == 0) {
				temp = CBSB_Sim(input[i], start, x0_pre[i]);
				vt_FG = get<1>(temp);
				FGInput_Steam[i] = { vt_FG[0] * (vt_FG[9] / 3.6) * (1.0 - vt_FG[10]),vt_FG[17] };
			}
			else {
				temp = CBSB_Simple_Sim(input[i]);
			}
		}
		else if (heater[i] == "CFB_Steam_Boiler") {
			if (*(input[i].end() - 1) == 0) {
				temp = CFBSB_Sim(input[i], start, x0_pre[i]);
				vt_FG = get<1>(temp);
				FGInput_Steam[i] = { vt_FG[0] * (vt_FG[9] / 3.6) * (1.0 - vt_FG[10]),vt_FG[17] };
			}
			else {
				temp = CFBSB_Simple_Sim(input[i]);
			}
		}
		else if (heater[i] == "Gas_Steam_Boiler") {
			if (*(input[i].end() - 1) == 0) {
				temp = GSB_Sim(input[i], start, x0_pre[i]);
				vt_FG = get<1>(temp);
				FGInput_Steam[i] = { vt_FG[0] * (vt_FG[9] / 3600.0) * (1.0 - vt_FG[10]),vt_FG[19] };
			}
			else {
				temp = GSB_Simple_Sim(input[i]);
			}
		}
		st = Add(get<0>(temp), "_heater_" + to_string(i + 1));
		vt = get<1>(temp);

		if (heater[i] == "Elec_Steam_Boiler") {
			x0_pre_i.push_back(vt[5]);
		}
		else if (heater[i] == "CB_Steam_Boiler" || heater[i] == "CFB_Steam_Boiler") {
			if (*(input[i].end() - 1) == 0) {
				x0_pre_i.push_back(vt[16]);
			}
			else {
				x0_pre_i.push_back(vt[9]);
			}
		}
		else if (heater[i] == "Gas_Steam_Boiler") {
			if (*(input[i].end() - 1) == 0) {
				x0_pre_i.push_back(vt[18]);
			}
			else {
				x0_pre_i.push_back(vt[9]);
			}
		}
		para_name.insert(para_name.end(), st.begin(), st.end());
		value.insert(value.end(), vt.begin(), vt.end());
		if (LB_St_Con.find(i) == LB_St_Con.end()) {
			if (heater[i] == "Elec_Steam_Boiler") {
				T_in[i] = vt[5];
				M_in[i] = vt[3] / 3.6; // kg/s
				P_in[i] = vt[6];
				H_in[i] = vt[8];
			}
			else if (heater[i] == "CB_Steam_Boiler" || heater[i] == "CFB_Steam_Boiler") {
				if (*(input[i].end() - 1) == 0) {
					T_in[i] = vt[15];
					M_in[i] = vt[2] / 3.6; // kg/s
					P_in[i] = vt[13];
					H_in[i] = vt[11];
				}
				else {
					T_in[i] = vt[9];
					M_in[i] = vt[4] / 3.6; // kg/s
					P_in[i] = vt[8];
					H_in[i] = vt[6];
				}
			}
			else if (heater[i] == "Gas_Steam_Boiler") {
				if (*(input[i].end() - 1) == 0) {
					T_in[i] = vt[17];
					M_in[i] = vt[2] / 3.6; // kg/s
					P_in[i] = vt[15];
					H_in[i] = vt[11];
				}
				else {
					T_in[i] = vt[9];
					M_in[i] = vt[4] / 3.6; // kg/s
					P_in[i] = vt[8];
					H_in[i] = vt[6];
				}
			}
		}
		else {
			T_in[i] = 0.0;
			M_in[i] = 0.0;
			P_in[i] = 0.0;
			H_in[i] = 0.0;

			if (heater[i] == "Elec_Steam_Boiler") {
				LB_Steam_part[LB_St_Con[i]].push_back({ vt[5], vt[3] / 3.6, vt[6], vt[8] });
			}
			else if (heater[i] == "CB_Steam_Boiler" || heater[i] == "CFB_Steam_Boiler") {
				if (*(input[i].end() - 1) == 0) {
					LB_Steam_part[LB_St_Con[i]].push_back({ vt[15], vt[2] / 3.6, vt[13], vt[11] });
				}
				else {
					LB_Steam_part[LB_St_Con[i]].push_back({ vt[9], vt[4] / 3.6, vt[8], vt[6] });
				}
			}
			else if (heater[i] == "Gas_Steam_Boiler") {
				if (*(input[i].end() - 1) == 0) {
					LB_Steam_part[LB_St_Con[i]].push_back({ vt[17], vt[2] / 3.6, vt[15], vt[11] });
				}
				else {
					LB_Steam_part[LB_St_Con[i]].push_back({ vt[9], vt[4] / 3.6, vt[8], vt[6] });
				}
			}
		}
		x0_pre_all.push_back(x0_pre_i);
	}

	for (int is = 0; is < LB_num; ++is) {
		int len_local = LB_Steam_part[is].size();
		vector<double> T_in_temp(len_local);
		vector<double> H_in_temp(len_local);
		vector<double> M_in_temp(len_local);
		vector<double> P_in_temp(len_local);
		for (int iss = 0; iss < len_local; ++iss) {
			T_in_temp[iss] = LB_Steam_part[is][iss][0];
			H_in_temp[iss] = LB_Steam_part[is][iss][3];
			M_in_temp[iss] = LB_Steam_part[is][iss][1];
			P_in_temp[iss] = LB_Steam_part[is][iss][2];
		}
		temp = SteamHub_Sim(T_in_temp, H_in_temp, M_in_temp, P_in_temp);
		Steam2HVAC[is] = { *(get<1>(temp).end() - 4), *(get<1>(temp).end() - 3) * 3.6 };
	}
	if (!start) {
		for (int i = 0; i < FGWasteHeatUse_S[0].size() + FGWasteHeatUse_S[1].size(); ++i) {
			vector<double> input_WHSB_use;
			input_WHSB_use.assign(input_WHSB[i].begin(), input_WHSB[i].end());
			if (i < FGWasteHeatUse_S[0].size()) {
				input_WHSB_use[0] = FGInput_Water[FGWasteHeatUse_S[0][i]][0];
			}
			else {
				input_WHSB_use[0] = FGInput_Steam[FGWasteHeatUse_S[1][i - FGWasteHeatUse_S[0].size()]][0];
			}
			temp_WHB = WHSB_Sim(input_WHSB_use);

			vs_WHB_temp.assign(get<0>(temp_WHB).begin(), get<0>(temp_WHB).end());
			vs_WHB_temp = Add(vs_WHB_temp, "_WasteHeatBoiler_" + to_string(i + 1));

			vs_WHB.insert(vs_WHB.end(), vs_WHB_temp.begin(), vs_WHB_temp.end());
			vt_WHB.insert(vt_WHB.end(), get<1>(temp_WHB).begin(), get<1>(temp_WHB).end());

			T_in.push_back(get<1>(temp_WHB)[6]); // K
			H_in.push_back(get<1>(temp_WHB)[4]); // kJ/kg
			M_in.push_back(get<1>(temp_WHB)[2]); // kg/s
			P_in.push_back(get<1>(temp_WHB)[5]); // MPa
		}
	}
	temp = SteamHub_Sim(T_in, H_in, M_in, P_in);
	st = get<0>(temp);
	vt = get<1>(temp);

	para_name.insert(para_name.end(), vs_WHB.begin(), vs_WHB.end());
	value.insert(value.end(), vt_WHB.begin(), vt_WHB.end());

	para_name.insert(para_name.end(), st.begin(), st.end());
	value.insert(value.end(), vt.begin(), vt.end());

	return{ para_name,value,x0_pre_all,Steam2HVAC };
}

tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> CWSystemSim_single_step_steam(double dt, const vector<string>& chiller, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, const vector<double>& T0, vector<int> charge, bool discharge, bool start, vector<vector<double>> x0_pre, vector<vector<double>>& Steam_in_out, vector<string> LB_TP) {

	if (!with_storage || discharge) { fill(charge.begin(), charge.end(), 0.0); }

	vector<string> para_name;
	vector<double> value;
	int len = chiller.size();
	vector<double> T_in(len), M_in(len);

	tuple<vector<string>, vector<double>> temp, temp_FGHub;
	tuple<vector<string>, vector<double>, vector<double>> temp_ST;
	vector<string> st, vs_temp_FG;
	vector<double> vt, vv_temp_FG;
	vector<double> T_in_ST, M_in_ST;
	vector<vector<double>> x0_pre_all;

	for (int i = 0; i < len; ++i) {
		vector<double> x0_pre_i;

		if (chiller[i] == "ACC") {
			temp = ACC_Sim(input[i], start, x0_pre[i]);
		}
		else if (chiller[i] == "WCC") {
			temp = WCC_Sim(input[i], start, x0_pre[i]);
		}
		else if (chiller[i] == "GSHP") {
			temp = GSHP_Sim(input[i], "Cooling", start, x0_pre[i]);
		}
		else if (chiller[i] == "LiBr_Steam") {
			temp = LiBr_Steam_Sim(input[i], start, x0_pre[i], Steam_in_out[i], "Cooling", LB_TP[i]);
		}
		else if (chiller[i] == "LiBr_DC") {
			temp = LiBr_DC_Sim(input[i], start, x0_pre[i], "Cooling");
		}
		else if (chiller[i] == "LiBr_WHW") {
			temp = LiBr_WHW_Sim(input[i], start, x0_pre[i], "Cooling");
		}
		else if (chiller[i] == "LiBr_Solar") {
			temp = LiBr_Solar_Sim(input[i], start, x0_pre[i], "Cooling", dt);
		}
		else if (chiller[i] == "LiBr_FG") {
			vector<double> Qin_FG_temp, Tin_FG_temp;
			int current_index_LBFG = LB_FG_INDEX.front();
			string current_type_LBFG = LB_FG_TYPE.front();
			if (get<0>(LB_FG_Source[current_index_LBFG]) == "Steam") {
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					Qin_FG_temp.push_back(FGInput_Steam[ii][0]);
					Tin_FG_temp.push_back(FGInput_Steam[ii][1]);
				}
			}
			else if (get<0>(LB_FG_Source[current_index_LBFG]) == "Water") {
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					Qin_FG_temp.push_back(FGInput_Water[ii][0]);
					Tin_FG_temp.push_back(FGInput_Water[ii][1]);
				}
			}
			else if (get<0>(LB_FG_Source[current_index_LBFG]) == "Steam+Water") {
				int countt = -1;
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					countt++;
					if (ii == -1) { break; }
					else {
						Qin_FG_temp.push_back(FGInput_Steam[ii][0]);
						Tin_FG_temp.push_back(FGInput_Steam[ii][1]);
					}
				}
				vector<int> WB_QT;
				WB_QT.assign(get<1>(LB_FG_Source[current_index_LBFG]).begin() + countt + 1, get<1>(LB_FG_Source[current_index_LBFG]).end());
				for (int ii : WB_QT) {
					Qin_FG_temp.push_back(FGInput_Water[ii][0]);
					Tin_FG_temp.push_back(FGInput_Water[ii][1]);
				}
			}
			temp_FGHub = FGHub_Sim(Qin_FG_temp, Tin_FG_temp);
			vector<double> input_FG;
			input_FG.assign(get<1>(temp_FGHub).end() - 4, get<1>(temp_FGHub).end());
			temp = LiBr_FG_Sim(input[i], start, x0_pre[i], "Cooling", current_type_LBFG, input_FG);
			vs_temp_FG = get<0>(temp);
			vv_temp_FG = get<1>(temp);
			vs_temp_FG.insert(vs_temp_FG.end(), get<0>(temp_FGHub).begin(), get<0>(temp_FGHub).end());
			vv_temp_FG.insert(vv_temp_FG.end(), get<1>(temp_FGHub).begin(), get<1>(temp_FGHub).end());
			temp = { vs_temp_FG,vv_temp_FG };

			LB_FG_INDEX.pop_front();
			LB_FG_TYPE.pop_front();
			LB_FG_INDEX.push_back(current_index_LBFG);
			LB_FG_TYPE.push_back(current_type_LBFG);
		}
		st = Add(get<0>(temp), "_chiller_" + to_string(i + 1));
		vt = get<1>(temp);

		x0_pre_i.push_back(vt[17]);

		para_name.insert(para_name.end(), st.begin(), st.end());
		value.insert(value.end(), vt.begin(), vt.end());

		T_in[i] = vt[3] * (1.0 - charge[i]);
		M_in[i] = vt[21] * (1.0 - charge[i]);

		if (charge[i]) {
			T_in_ST.push_back(vt[3]);
			M_in_ST.push_back(vt[21]);
		}

		x0_pre_all.push_back(x0_pre_i);
	}

	if (with_storage)
	{
		if (discharge) {
			temp_ST = STank_Sim(input_ST, T0, dt, "Cooling");
			T_in.push_back(get<1>(temp_ST)[9]);
			M_in.push_back(-get<1>(temp_ST)[0]);
		}
		else {
			vector<double> input_ST_use = input_ST;
			if (!T_in_ST.empty()) {
				vector<double> input_ST_part = get<1>(WHub_Sim(T_in_ST, M_in_ST));
				input_ST_use[0] = input_ST_part[input_ST_part.size() - 1];
				input_ST_use[1] = input_ST_part[input_ST_part.size() - 2];
			}
			else {
				input_ST_use[0] = 0.0;
				input_ST_use[1] = 0.0;
			}

			temp_ST = STank_Sim(input_ST_use, T0, dt, "Cooling");
			T_in.push_back(0.0);
			M_in.push_back(0.0);
		}
		st = get<0>(temp_ST);
		vt = get<1>(temp_ST);
		para_name.insert(para_name.end(), st.begin(), st.end());
		value.insert(value.end(), vt.begin(), vt.end());
	}

	temp = WHub_Sim(T_in, M_in);
	st = get<0>(temp);
	vt = get<1>(temp);
	para_name.insert(para_name.end(), st.begin(), st.end());
	value.insert(value.end(), vt.begin(), vt.end());

	if (with_storage) { return{ para_name,value,get<2>(temp_ST),x0_pre_all }; }
	else { return{ para_name,value,{0},x0_pre_all }; }
}

tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> CWSystemSim_single_step(double dt, const vector<string>& chiller, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, const vector<double>& T0, vector<int> charge, bool discharge, bool start, vector<vector<double>> x0_pre) {

	if (!with_storage || discharge) { fill(charge.begin(), charge.end(), 0.0); }

	vector<string> para_name;
	vector<double> value;
	int len = chiller.size();
	vector<double> T_in(len), M_in(len);

	tuple<vector<string>, vector<double>> temp, temp_FGHub;
	tuple<vector<string>, vector<double>, vector<double>> temp_ST;
	vector<string> st, vs_temp_FG;
	vector<double> vt, vv_temp_FG;
	vector<double> T_in_ST, M_in_ST;
	vector<vector<double>> x0_pre_all;

	for (int i = 0; i < len; ++i) {
		vector<double> x0_pre_i;

		if (chiller[i] == "ACC") {
			temp = ACC_Sim(input[i], start, x0_pre[i]);
		}
		else if (chiller[i] == "WCC") {
			temp = WCC_Sim(input[i], start, x0_pre[i]);
		}
		else if (chiller[i] == "GSHP") {
			temp = GSHP_Sim(input[i], "Cooling", start, x0_pre[i]);
		}
		else if (chiller[i] == "LiBr_DC") {
			temp = LiBr_DC_Sim(input[i], start, x0_pre[i], "Cooling");
		}
		else if (chiller[i] == "LiBr_WHW") {
			temp = LiBr_WHW_Sim(input[i], start, x0_pre[i], "Cooling");
		}
		else if (chiller[i] == "LiBr_Solar") {
			temp = LiBr_Solar_Sim(input[i], start, x0_pre[i], "Cooling", dt);
		}
		else if (chiller[i] == "LiBr_FG") {
			vector<double> Qin_FG_temp, Tin_FG_temp;
			int current_index_LBFG = LB_FG_INDEX.front();
			string current_type_LBFG = LB_FG_TYPE.front();
			if (get<0>(LB_FG_Source[current_index_LBFG]) == "Steam") {
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					Qin_FG_temp.push_back(FGInput_Steam[ii][0]);
					Tin_FG_temp.push_back(FGInput_Steam[ii][1]);
				}
			}
			else if (get<0>(LB_FG_Source[current_index_LBFG]) == "Water") {
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					Qin_FG_temp.push_back(FGInput_Water[ii][0]);
					Tin_FG_temp.push_back(FGInput_Water[ii][1]);
				}
			}
			else if (get<0>(LB_FG_Source[current_index_LBFG]) == "Steam+Water") {
				int countt = -1;
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					countt++;
					if (ii == -1) { break; }
					else {
						Qin_FG_temp.push_back(FGInput_Steam[ii][0]);
						Tin_FG_temp.push_back(FGInput_Steam[ii][1]);
					}
				}
				vector<int> WB_QT;
				WB_QT.assign(get<1>(LB_FG_Source[current_index_LBFG]).begin() + countt + 1, get<1>(LB_FG_Source[current_index_LBFG]).end());
				for (int ii : WB_QT) {
					Qin_FG_temp.push_back(FGInput_Water[ii][0]);
					Tin_FG_temp.push_back(FGInput_Water[ii][1]);
				}
			}
			temp_FGHub = FGHub_Sim(Qin_FG_temp, Tin_FG_temp);
			vector<double> input_FG;
			input_FG.assign(get<1>(temp_FGHub).end() - 4, get<1>(temp_FGHub).end());
			temp = LiBr_FG_Sim(input[i], start, x0_pre[i], "Cooling", current_type_LBFG, input_FG);
			vs_temp_FG = get<0>(temp);
			vv_temp_FG = get<1>(temp);
			vs_temp_FG.insert(vs_temp_FG.end(), get<0>(temp_FGHub).begin(), get<0>(temp_FGHub).end());
			vv_temp_FG.insert(vv_temp_FG.end(), get<1>(temp_FGHub).begin(), get<1>(temp_FGHub).end());
			temp = { vs_temp_FG,vv_temp_FG };

			LB_FG_INDEX.pop_front();
			LB_FG_TYPE.pop_front();
			LB_FG_INDEX.push_back(current_index_LBFG);
			LB_FG_TYPE.push_back(current_type_LBFG);
		}
		st = Add(get<0>(temp), "_chiller_" + to_string(i + 1));
		vt = get<1>(temp);

		x0_pre_i.push_back(vt[17]);

		para_name.insert(para_name.end(), st.begin(), st.end());
		value.insert(value.end(), vt.begin(), vt.end());

		T_in[i] = vt[3] * (1.0 - charge[i]);
		M_in[i] = vt[21] * (1.0 - charge[i]);

		if (charge[i]) {
			T_in_ST.push_back(vt[3]);
			M_in_ST.push_back(vt[21]);
		}

		x0_pre_all.push_back(x0_pre_i);
	}

	if (with_storage)
	{
		if (discharge) {
			temp_ST = STank_Sim(input_ST, T0, dt, "Cooling");
			T_in.push_back(get<1>(temp_ST)[9]);
			M_in.push_back(-get<1>(temp_ST)[0]);
		}
		else {
			vector<double> input_ST_use = input_ST;
			if (!T_in_ST.empty()) {
				vector<double> input_ST_part = get<1>(WHub_Sim(T_in_ST, M_in_ST));
				input_ST_use[0] = input_ST_part[input_ST_part.size() - 1];
				input_ST_use[1] = input_ST_part[input_ST_part.size() - 2];
			}
			else {
				input_ST_use[0] = 0.0;
				input_ST_use[1] = 0.0;
			}

			temp_ST = STank_Sim(input_ST_use, T0, dt, "Cooling");
			T_in.push_back(0.0);
			M_in.push_back(0.0);
		}
		st = get<0>(temp_ST);
		vt = get<1>(temp_ST);
		para_name.insert(para_name.end(), st.begin(), st.end());
		value.insert(value.end(), vt.begin(), vt.end());
	}

	temp = WHub_Sim(T_in, M_in);
	st = get<0>(temp);
	vt = get<1>(temp);
	para_name.insert(para_name.end(), st.begin(), st.end());
	value.insert(value.end(), vt.begin(), vt.end());

	if (with_storage) { return{ para_name,value,get<2>(temp_ST),x0_pre_all }; }
	else { return{ para_name,value,{0},x0_pre_all }; }
}

tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> HWSystemSim_single_step_steam(double dt, const vector<string>& heater, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, const vector<double>& T0, vector<int> charge, bool discharge, bool start, vector<vector<double>> x0_pre, vector<vector<double>>& Steam_in_out, vector<string> LB_TP, vector<vector<double>> FGWasteHeatUse, vector<vector<double>> input_WHHWB) {

	if (!with_storage || discharge) { fill(charge.begin(), charge.end(), 0.0); }

	vector<string> para_name;
	vector<double> value;
	int len = heater.size();
	vector<double> T_in(len), M_in(len);

	tuple<vector<string>, vector<double>> temp, temp_FGHub, temp_WHB;
	tuple<vector<string>, vector<double>, vector<double>> temp_ST;
	vector<string> st, vs_temp_FG, vs_WHB, vs_WHB_temp;
	vector<double> vt, vt_FG, vv_temp_FG, vt_WHB;
	vector<double> T_in_ST, M_in_ST;
	vector<vector<double>> x0_pre_all;

	for (int i = 0; i < len; ++i) {
		vector<double> x0_pre_i;

		if (heater[i] == "ASHP") {
			temp = ASHP_Sim(input[i], start, x0_pre[i]);
		}
		else if (heater[i] == "WSHP") {
			temp = WSHP_Sim(input[i], start, x0_pre[i]);
		}
		else if (heater[i] == "GSHP") {
			temp = GSHP_Sim(input[i], "Heating", start, x0_pre[i]);
		}
		else if (heater[i] == "LiBr_Steam") {
			temp = LiBr_Steam_Sim(input[i], start, x0_pre[i], Steam_in_out[i], "Heating", LB_TP[i]);
		}
		else if (heater[i] == "LiBr_DC") {
			temp = LiBr_DC_Sim(input[i], start, x0_pre[i], "Heating");
		}
		else if (heater[i] == "LiBr_WHW") {
			temp = LiBr_WHW_Sim(input[i], start, x0_pre[i], "Heating");
		}
		else if (heater[i] == "LiBr_Solar") {
			temp = LiBr_Solar_Sim(input[i], start, x0_pre[i], "Heating", dt);
		}
		else if (heater[i] == "CB_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				temp = CBHWB_Sim(input[i], start, x0_pre[i]);
				vt_FG = get<1>(temp);
				FGInput_Water[i] = { vt_FG[0] * (vt_FG[9] / 3.6) * (1.0 - vt_FG[10]),vt_FG[13] };
			}
			else {
				temp = CBHWB_Simple_Sim(input[i]);
			}
		}
		else if (heater[i] == "CFB_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				temp = CFBHWB_Sim(input[i], start, x0_pre[i]);
				vt_FG = get<1>(temp);
				FGInput_Water[i] = { vt_FG[0] * (vt_FG[9] / 3.6) * (1.0 - vt_FG[10]),vt_FG[13] };
			}
			else {
				temp = CFBHWB_Simple_Sim(input[i]);
			}
		}
		else if (heater[i] == "Gas_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				temp = GHWB_Sim(input[i], start, x0_pre[i]);
				vt_FG = get<1>(temp);
				FGInput_Water[i] = { vt_FG[0] * (vt_FG[9] / 3600.0) * (1.0 - vt_FG[10]),vt_FG[15] };
			}
			else {
				temp = GHWB_Simple_Sim(input[i]);
			}
		}
		else if (heater[i] == "Elec_Boiler") {
			temp = EHWB_Sim(input[i]);
		}
		else if (heater[i] == "LiBr_FG") {
			vector<double> Qin_FG_temp, Tin_FG_temp;
			int current_index_LBFG = LB_FG_INDEX.front();
			string current_type_LBFG = LB_FG_TYPE.front();
			if (get<0>(LB_FG_Source[current_index_LBFG]) == "Steam") {
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					Qin_FG_temp.push_back(FGInput_Steam[ii][0]);
					Tin_FG_temp.push_back(FGInput_Steam[ii][1]);
				}
			}
			else if (get<0>(LB_FG_Source[current_index_LBFG]) == "Water") {
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					Qin_FG_temp.push_back(FGInput_Water[ii][0]);
					Tin_FG_temp.push_back(FGInput_Water[ii][1]);
				}
			}
			else if (get<0>(LB_FG_Source[current_index_LBFG]) == "Steam+Water") {
				int countt = -1;
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					countt++;
					if (ii == -1) { break; }
					else {
						Qin_FG_temp.push_back(FGInput_Steam[ii][0]);
						Tin_FG_temp.push_back(FGInput_Steam[ii][1]);
					}
				}
				vector<int> WB_QT;
				WB_QT.assign(get<1>(LB_FG_Source[current_index_LBFG]).begin() + countt + 1, get<1>(LB_FG_Source[current_index_LBFG]).end());
				for (int ii : WB_QT) {
					Qin_FG_temp.push_back(FGInput_Water[ii][0]);
					Tin_FG_temp.push_back(FGInput_Water[ii][1]);
				}
			}
			temp_FGHub = FGHub_Sim(Qin_FG_temp, Tin_FG_temp);
			vector<double> input_FG;
			input_FG.assign(get<1>(temp_FGHub).end() - 4, get<1>(temp_FGHub).end());
			temp = LiBr_FG_Sim(input[i], start, x0_pre[i], "Heating", current_type_LBFG, input_FG);
			vs_temp_FG = get<0>(temp);
			vv_temp_FG = get<1>(temp);
			vs_temp_FG.insert(vs_temp_FG.end(), get<0>(temp_FGHub).begin(), get<0>(temp_FGHub).end());
			vv_temp_FG.insert(vv_temp_FG.end(), get<1>(temp_FGHub).begin(), get<1>(temp_FGHub).end());
			temp = { vs_temp_FG,vv_temp_FG };

			LB_FG_INDEX.pop_front();
			LB_FG_TYPE.pop_front();
			LB_FG_INDEX.push_back(current_index_LBFG);
			LB_FG_TYPE.push_back(current_type_LBFG);
		}
		st = Add(get<0>(temp), "_heater_" + to_string(i + 1));
		vt = get<1>(temp);

		if (heater[i] == "CB_Boiler" || heater[i] == "CFB_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				x0_pre_i.push_back(vt[12]);
			}
			else {
				x0_pre_i.push_back(vt[6]);
			}
		}
		else if (heater[i] == "Gas_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				x0_pre_i.push_back(vt[14]);
			}
			else {
				x0_pre_i.push_back(vt[6]);
			}
		}
		else if (heater[i] == "Elec_Boiler") {
			x0_pre_i.push_back(vt[5]);
		}
		else {
			x0_pre_i.push_back(vt[17]);
		}

		para_name.insert(para_name.end(), st.begin(), st.end());
		value.insert(value.end(), vt.begin(), vt.end());
		if (heater[i] == "CB_Boiler" || heater[i] == "CFB_Boiler" || heater[i] == "Gas_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				T_in[i] = vt[11] * (1.0 - charge[i]);
				M_in[i] = vt[2] * (1.0 - charge[i]) / 3.6;
			}
			else {
				T_in[i] = vt[6] * (1.0 - charge[i]);
				M_in[i] = vt[4] * (1.0 - charge[i]) / 3.6;
			}
		}
		else if (heater[i] == "Elec_Boiler") {
			T_in[i] = vt[5] * (1.0 - charge[i]);
			M_in[i] = vt[3] * (1.0 - charge[i]) / 3.6;
		}
		else {
			T_in[i] = vt[3] * (1.0 - charge[i]);
			M_in[i] = vt[21] * (1.0 - charge[i]);
		}

		if (charge[i]) {
			if (heater[i] == "CB_Boiler" || heater[i] == "CFB_Boiler" || heater[i] == "Gas_Boiler") {
				if (*(input[i].end() - 1) == 0.0) {
					T_in_ST.push_back(vt[11]);
					M_in_ST.push_back(vt[2] / 3.6);
				}
				else {
					T_in_ST.push_back(vt[6]);
					M_in_ST.push_back(vt[4] / 3.6);
				}
			}
			else if (heater[i] == "Elec_Boiler") {
				T_in_ST.push_back(vt[5]);
				M_in_ST.push_back(vt[3] / 3.6);
			}
			else {
				T_in_ST.push_back(vt[3]);
				M_in_ST.push_back(vt[21]);
			}
		}

		x0_pre_all.push_back(x0_pre_i);
	}

	for (int i = 0; i < FGWasteHeatUse[0].size() + FGWasteHeatUse[1].size(); ++i) {
		vector<double> input_WHHWB_use;
		input_WHHWB_use.assign(input_WHHWB[i].begin(), input_WHHWB[i].end() - 1);
		if (i < FGWasteHeatUse[0].size()) {
			input_WHHWB_use[0] = FGInput_Water[FGWasteHeatUse[0][i]][0];
		}
		else {
			input_WHHWB_use[0] = FGInput_Steam[FGWasteHeatUse[1][i - FGWasteHeatUse[0].size()]][0];
		}
		temp_WHB = WHHWB_Sim(input_WHHWB_use);

		vs_WHB_temp.assign(get<0>(temp_WHB).begin(), get<0>(temp_WHB).end());
		vs_WHB_temp = Add(vs_WHB_temp, "_WasteHeatBoiler_" + to_string(i + 1));

		vs_WHB.insert(vs_WHB.end(), vs_WHB_temp.begin(), vs_WHB_temp.end());
		vt_WHB.insert(vt_WHB.end(), get<1>(temp_WHB).begin(), get<1>(temp_WHB).end());

		if (*(input_WHHWB[i].end() - 1)) {
			T_in.push_back(0.0);
			M_in.push_back(0.0);
			T_in_ST.push_back(get<1>(temp_WHB)[4]);
			M_in_ST.push_back(get<1>(temp_WHB)[2]);
		}
		else {
			T_in.push_back(get<1>(temp_WHB)[4]);
			M_in.push_back(get<1>(temp_WHB)[2]);
		}
	}

	if (with_storage)
	{
		if (discharge) {
			temp_ST = STank_Sim(input_ST, T0, dt, "Heating");
			T_in.push_back(get<1>(temp_ST)[9]);
			M_in.push_back(get<1>(temp_ST)[0]);
		}
		else {
			vector<double> input_ST_use = input_ST;
			if (!T_in_ST.empty()) {
				vector<double> input_ST_part = get<1>(WHub_Sim(T_in_ST, M_in_ST));
				input_ST_use[0] = -input_ST_part[input_ST_part.size() - 1];
				input_ST_use[1] = input_ST_part[input_ST_part.size() - 2];
			}
			else {
				input_ST_use[0] = 0.0;
				input_ST_use[1] = 0.0;
			}

			temp_ST = STank_Sim(input_ST_use, T0, dt, "Heating");
			T_in.push_back(0.0);
			M_in.push_back(0.0);
		}
		st = get<0>(temp_ST);
		vt = get<1>(temp_ST);
		para_name.insert(para_name.end(), st.begin(), st.end());
		value.insert(value.end(), vt.begin(), vt.end());
	}

	temp = WHub_Sim(T_in, M_in);
	st = get<0>(temp);
	vt = get<1>(temp);

	para_name.insert(para_name.end(), vs_WHB.begin(), vs_WHB.end());
	value.insert(value.end(), vt_WHB.begin(), vt_WHB.end());

	para_name.insert(para_name.end(), st.begin(), st.end());
	value.insert(value.end(), vt.begin(), vt.end());

	if (with_storage) { return{ para_name,value,get<2>(temp_ST),x0_pre_all }; }
	else { return{ para_name,value,{0},x0_pre_all }; }
}

tuple<vector<string>, vector<double>, vector<double>, vector<vector<double>>> HWSystemSim_single_step(double dt, const vector<string>& heater, const vector<vector<double>>& input, bool with_storage, const vector<double>& input_ST, const vector<double>& T0, vector<int> charge, bool discharge, bool start, vector<vector<double>> x0_pre, vector<vector<double>> FGWasteHeatUse, vector<vector<double>> input_WHHWB) {

	if (!with_storage || discharge) { fill(charge.begin(), charge.end(), 0.0); }

	vector<string> para_name;
	vector<double> value;
	int len = heater.size();
	vector<double> T_in(len), M_in(len);

	tuple<vector<string>, vector<double>> temp, temp_FGHub, temp_WHB;
	tuple<vector<string>, vector<double>, vector<double>> temp_ST;
	vector<string> st, vs_temp_FG, vs_WHB, vs_WHB_temp;
	vector<double> vt, vt_FG, vv_temp_FG, vt_WHB;
	vector<double> T_in_ST, M_in_ST;
	vector<vector<double>> x0_pre_all;

	for (int i = 0; i < len; ++i) {
		vector<double> x0_pre_i;

		if (heater[i] == "ASHP") {
			temp = ASHP_Sim(input[i], start, x0_pre[i]);
		}
		else if (heater[i] == "WSHP") {
			temp = WSHP_Sim(input[i], start, x0_pre[i]);
		}
		else if (heater[i] == "GSHP") {
			temp = GSHP_Sim(input[i], "Heating", start, x0_pre[i]);
		}
		else if (heater[i] == "LiBr_DC") {
			temp = LiBr_DC_Sim(input[i], start, x0_pre[i], "Heating");
		}
		else if (heater[i] == "LiBr_WHW") {
			temp = LiBr_WHW_Sim(input[i], start, x0_pre[i], "Heating");
		}
		else if (heater[i] == "LiBr_Solar") {
			temp = LiBr_Solar_Sim(input[i], start, x0_pre[i], "Heating", dt);
		}
		else if (heater[i] == "CB_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				temp = CBHWB_Sim(input[i], start, x0_pre[i]);
				vt_FG = get<1>(temp);
				FGInput_Water[i] = { vt_FG[0] * (vt_FG[9] / 3.6) * (1.0 - vt_FG[10]),vt_FG[13] };
			}
			else {
				temp = CBHWB_Simple_Sim(input[i]);
			}
		}
		else if (heater[i] == "CFB_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				temp = CFBHWB_Sim(input[i], start, x0_pre[i]);
				vt_FG = get<1>(temp);
				FGInput_Water[i] = { vt_FG[0] * (vt_FG[9] / 3.6) * (1.0 - vt_FG[10]),vt_FG[13] };
			}
			else {
				temp = CFBHWB_Simple_Sim(input[i]);
			}
		}
		else if (heater[i] == "Gas_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				temp = GHWB_Sim(input[i], start, x0_pre[i]);
				vt_FG = get<1>(temp);
				FGInput_Water[i] = { vt_FG[0] * (vt_FG[9] / 3600.0) * (1.0 - vt_FG[10]),vt_FG[15] };
			}
			else {
				temp = GHWB_Simple_Sim(input[i]);
			}
		}
		else if (heater[i] == "Elec_Boiler") {
			temp = EHWB_Sim(input[i]);
		}
		else if (heater[i] == "LiBr_FG") {
			vector<double> Qin_FG_temp, Tin_FG_temp;
			int current_index_LBFG = LB_FG_INDEX.front();
			string current_type_LBFG = LB_FG_TYPE.front();
			if (get<0>(LB_FG_Source[current_index_LBFG]) == "Steam") {
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					Qin_FG_temp.push_back(FGInput_Steam[ii][0]);
					Tin_FG_temp.push_back(FGInput_Steam[ii][1]);
				}
			}
			else if (get<0>(LB_FG_Source[current_index_LBFG]) == "Water") {
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					Qin_FG_temp.push_back(FGInput_Water[ii][0]);
					Tin_FG_temp.push_back(FGInput_Water[ii][1]);
				}
			}
			else if (get<0>(LB_FG_Source[current_index_LBFG]) == "Steam+Water") {
				int countt = -1;
				for (int ii : get<1>(LB_FG_Source[current_index_LBFG])) {
					countt++;
					if (ii == -1) { break; }
					else {
						Qin_FG_temp.push_back(FGInput_Steam[ii][0]);
						Tin_FG_temp.push_back(FGInput_Steam[ii][1]);
					}
				}
				vector<int> WB_QT;
				WB_QT.assign(get<1>(LB_FG_Source[current_index_LBFG]).begin() + countt + 1, get<1>(LB_FG_Source[current_index_LBFG]).end());
				for (int ii : WB_QT) {
					Qin_FG_temp.push_back(FGInput_Water[ii][0]);
					Tin_FG_temp.push_back(FGInput_Water[ii][1]);
				}
			}
			temp_FGHub = FGHub_Sim(Qin_FG_temp, Tin_FG_temp);
			vector<double> input_FG;
			input_FG.assign(get<1>(temp_FGHub).end() - 4, get<1>(temp_FGHub).end());
			temp = LiBr_FG_Sim(input[i], start, x0_pre[i], "Heating", current_type_LBFG, input_FG);
			vs_temp_FG = get<0>(temp);
			vv_temp_FG = get<1>(temp);
			vs_temp_FG.insert(vs_temp_FG.end(), get<0>(temp_FGHub).begin(), get<0>(temp_FGHub).end());
			vv_temp_FG.insert(vv_temp_FG.end(), get<1>(temp_FGHub).begin(), get<1>(temp_FGHub).end());
			temp = { vs_temp_FG,vv_temp_FG };

			LB_FG_INDEX.pop_front();
			LB_FG_TYPE.pop_front();
			LB_FG_INDEX.push_back(current_index_LBFG);
			LB_FG_TYPE.push_back(current_type_LBFG);
		}
		st = Add(get<0>(temp), "_heater_" + to_string(i + 1));
		vt = get<1>(temp);

		if (heater[i] == "CB_Boiler" || heater[i] == "CFB_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				x0_pre_i.push_back(vt[12]);
			}
			else {
				x0_pre_i.push_back(vt[6]);
			}
		}
		else if (heater[i] == "Gas_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				x0_pre_i.push_back(vt[14]);
			}
			else {
				x0_pre_i.push_back(vt[6]);
			}
		}
		else if (heater[i] == "Elec_Boiler") {
			x0_pre_i.push_back(vt[5]);
		}
		else {
			x0_pre_i.push_back(vt[17]);
		}

		para_name.insert(para_name.end(), st.begin(), st.end());
		value.insert(value.end(), vt.begin(), vt.end());
		if (heater[i] == "CB_Boiler" || heater[i] == "CFB_Boiler" || heater[i] == "Gas_Boiler") {
			if (*(input[i].end() - 1) == 0.0) {
				T_in[i] = vt[11] * (1.0 - charge[i]);
				M_in[i] = vt[2] * (1.0 - charge[i]) / 3.6;
			}
			else {
				T_in[i] = vt[6] * (1.0 - charge[i]);
				M_in[i] = vt[4] * (1.0 - charge[i]) / 3.6;
			}
		}
		else if (heater[i] == "Elec_Boiler") {
			T_in[i] = vt[5] * (1.0 - charge[i]);
			M_in[i] = vt[3] * (1.0 - charge[i]) / 3.6;
		}
		else {
			T_in[i] = vt[3] * (1.0 - charge[i]);
			M_in[i] = vt[21] * (1.0 - charge[i]);
		}

		if (charge[i]) {
			if (heater[i] == "CB_Boiler" || heater[i] == "CFB_Boiler" || heater[i] == "Gas_Boiler") {
				if (*(input[i].end() - 1) == 0.0) {
					T_in_ST.push_back(vt[11]);
					M_in_ST.push_back(vt[2] / 3.6);
				}
				else {
					T_in_ST.push_back(vt[6]);
					M_in_ST.push_back(vt[4] / 3.6);
				}
			}
			else if (heater[i] == "Elec_Boiler") {
				T_in_ST.push_back(vt[5]);
				M_in_ST.push_back(vt[3] / 3.6);
			}
			else {
				T_in_ST.push_back(vt[3]);
				M_in_ST.push_back(vt[21]);
			}
		}

		x0_pre_all.push_back(x0_pre_i);
	}

	for (int i = 0; i < FGWasteHeatUse[0].size(); ++i) {
		vector<double> input_WHHWB_use;
		input_WHHWB_use.assign(input_WHHWB[i].begin(), input_WHHWB[i].end() - 1);
		input_WHHWB_use[0] = FGInput_Water[FGWasteHeatUse[0][i]][0];
		temp_WHB = WHHWB_Sim(input_WHHWB_use);
		
		vs_WHB_temp.assign(get<0>(temp_WHB).begin(), get<0>(temp_WHB).end());
		vs_WHB_temp = Add(vs_WHB_temp, "_WasteHeatBoiler_" + to_string(i + 1));
		
		vs_WHB.insert(vs_WHB.end(), vs_WHB_temp.begin(), vs_WHB_temp.end());
		vt_WHB.insert(vt_WHB.end(), get<1>(temp_WHB).begin(), get<1>(temp_WHB).end());
		
		if (*(input_WHHWB[i].end() - 1)) {
			T_in.push_back(0.0);
			M_in.push_back(0.0);
			T_in_ST.push_back(get<1>(temp_WHB)[4]);
			M_in_ST.push_back(get<1>(temp_WHB)[2]);
		}
		else {
			T_in.push_back(get<1>(temp_WHB)[4]);
			M_in.push_back(get<1>(temp_WHB)[2]);
		}
	}

	if (with_storage)
	{
		if (discharge) {
			temp_ST = STank_Sim(input_ST, T0, dt, "Heating");
			T_in.push_back(get<1>(temp_ST)[9]);
			M_in.push_back(get<1>(temp_ST)[0]);
		}
		else {
			vector<double> input_ST_use = input_ST;
			if (!T_in_ST.empty()) {
				vector<double> input_ST_part = get<1>(WHub_Sim(T_in_ST, M_in_ST));
				input_ST_use[0] = -input_ST_part[input_ST_part.size() - 1];
				input_ST_use[1] = input_ST_part[input_ST_part.size() - 2];
			}
			else {
				input_ST_use[0] = 0.0;
				input_ST_use[1] = 0.0;
			}

			temp_ST = STank_Sim(input_ST_use, T0, dt, "Heating");
			T_in.push_back(0.0);
			M_in.push_back(0.0);
		}
		st = get<0>(temp_ST);
		vt = get<1>(temp_ST);
		para_name.insert(para_name.end(), st.begin(), st.end());
		value.insert(value.end(), vt.begin(), vt.end());
	}

	temp = WHub_Sim(T_in, M_in);
	st = get<0>(temp);
	vt = get<1>(temp);
	
	para_name.insert(para_name.end(), vs_WHB.begin(), vs_WHB.end());
	value.insert(value.end(), vt_WHB.begin(), vt_WHB.end());

	para_name.insert(para_name.end(), st.begin(), st.end());
	value.insert(value.end(), vt.begin(), vt.end());

	if (with_storage) { return{ para_name,value,get<2>(temp_ST),x0_pre_all }; }
	else { return{ para_name,value,{0},x0_pre_all }; }
}