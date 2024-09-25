#ifndef AE_Solve_Algorithm
#define AE_Solve_Algorithm
#include <vector>
#include <iostream>
#include <math.h>
#include <random>
#include <algorithm>
#include <functional>
using namespace std;

vector<vector<double>> AE_Solve_AHS(vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), double range_abs_max, vector<double> x0, vector<double> K, bool use_PSO);

vector<vector<double>> AE_Solve_Basic(vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), vector<double> x0, vector<double> K);

vector<vector<double>> AE_Solve_PSO(vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), double range_abs_max, int dim, double e, vector<double> K, double pos);

#endif