#ifndef ODE_Solve_Algorithm
#define ODE_Solve_Algorithm
#include <vector>
#include <iostream>
#include <math.h>
#include <random>
#include <algorithm>
#include <functional>
using namespace std;

vector<double> Scale(const vector<double>& v, double k);

vector<double> Add(const vector<double>& v1, const vector<double>& v2);

vector<string> Add(const vector<string>& v, string add);

vector<double> Substract(const vector<double>& v1, const vector<double>& v2);

vector<double> Square(const vector<double>& v);

vector<double> Multiply(const vector<double>& v1, const vector<double>& v2);

vector<double> Div(const vector<double>& v1, const vector<double>& v2);

double Max(const vector<double>& v);

double Tolerance(const vector<double>& v1, const vector<double>& v2);

vector<vector <double>> ODE_Solve_Euler_Forward(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_Euler_Trap(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_Euler_Backward(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_Adamas_Explicit(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K);

vector<vector <double>> ODE_Solve_Adamas_Implicit(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K);

vector<vector <double>> ODE_Solve_Adamas_PECE(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_Hamming(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_Nystrom_Explicit(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K);

vector<vector <double>> ODE_Solve_Milne_Simpson(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_MP2(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Heun2(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Heun3(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Kutta3(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Nystrom3(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Kutta4(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_C4(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKBS(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Gill4(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Nystrom5(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Luther5(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Butcher5(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Butcher5a(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Butcher5b(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Butcher5c(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKF45(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKF78(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Butcher6(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKDP(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKCK(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKHH(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK7S(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK6M(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_BDF(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K);

vector<vector <double>> ODE_Solve_NDF(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_Implicit(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K);

vector<vector <double>> ODE_Solve_Gear(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_RIA(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_RIIA(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_LIIIA(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_LIIIB(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_LIIIC(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_LIIICn(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_SDIRK(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K);

vector<vector <double>> ODE_Solve_RK_ESDIRK(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKK(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKV(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKT5(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

vector<vector <double>> ODE_Solve_RKT65(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K);

#endif