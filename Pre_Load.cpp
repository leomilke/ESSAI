#include "Sub_Models.h"
#include "ODE_Solve_Algorithm.h"
#include "AE_Solve_Algorithm.h"
#include "Pre_load.h"

vector<double> T_Steam_P(10000);
vector<vector<double>> h_Steam_Px(1000, vector<double>(100));
vector<vector<double>> s_Steam_Px(1000, vector<double>(100));
vector<vector<double>> T_Steam_hs(5000, vector<double>(1000));
bool read_steam = false;

void Read_Steam_Dat() {
    if (!read_steam) {
        ifstream file("Steam_T(P).dat");
        for (int i = 0; i < 10000; ++i) {
            file >> T_Steam_P[i];
        }
        file.close();

        file.open("Steam_h(P,x).dat");
        for (int i = 0; i < 1000; ++i) {
            for (int j = 0; j < 100; ++j) {
                file >> h_Steam_Px[i][j];
            }
        }
        file.close();

        file.open("Steam_s(P,x).dat");
        for (int i = 0; i < 1000; ++i) {
            for (int j = 0; j < 100; ++j) {
                file >> s_Steam_Px[i][j];
            }
        }
        file.close();

        file.open("Steam_T(h,s).dat");
        for (int i = 0; i < 5000; ++i) {
            for (int j = 0; j < 1000; ++j) {
                file >> T_Steam_hs[i][j];
            }
        }
        file.close();

        read_steam = true;
    }
}

double x_Ph(double P, double h) {
    double e_min = 1.0e10;
    double e, x;
    for (double i = 0.0; i < 100.0; ++i) {
        e = abs(h - h_Steam_Px[round(P * 100.0) - 1.0][i]);
        if (e < e_min) {
            e_min = e;
            x = (i + 1.0) / 100.0;
        }
    }

    return x;
}

double s_Px(double P, double x) {
    return s_Steam_Px[round(P * 100.0) - 1.0][round(x * 100.0) - 1.0];
}

double T_hs(double h, double s) {
    return T_Steam_hs[round(h) - 1.0][round(s * 100.0) - 1.0];
}

double P_xh(double x, double h) {
    double e_min = 1.0e10;
    double e, P;
    for (double i = 0.0; i < 1000.0; ++i) {
        e = abs(h - h_Steam_Px[i][round(x * 100.0) - 1.0]);
        if (e < e_min) {
            e_min = e;
            P = (i + 1.0) / 100.0;
        }
    }

    return P;
}

double Cal_T_Steam_P(double p) { // p(MPa)
    return T_Steam_P[round(p * 1000.0) - 1.0]; // unit: K
}

double Cal_P_Steam_T(double T) { // K
    double emin = 1.0e10, e, P;
    for (double p = 1; p <= 10000; ++p) {
        e = abs(T - T_Steam_P[p - 1.0]);
        if (e < emin) {
            emin = e;
            P = p / 1000.0;
        }
    }
    
    return P; // MPa
}

double Cal_T_FG(double h) { // kJ/kg
    return 0.96981 * h - 40.846; // T_fluegas (K)
}

double Cal_h_FG(double T) { // K
    return (T + 40.846) / 0.96981; // h, kJ/kg
}