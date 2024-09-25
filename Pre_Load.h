#ifndef Pre_Load
#define Pre_Load
#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

void Read_Steam_Dat();

double x_Ph(double P, double h);

double P_xh(double x, double h);

double T_hs(double h, double s);

double s_Px(double P, double x);

double Cal_T_Steam_P(double p);

double Cal_P_Steam_T(double T);

double Cal_T_FG(double h);

double Cal_h_FG(double T);

#endif