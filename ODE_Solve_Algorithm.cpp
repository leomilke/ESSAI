#include "ODE_Solve_Algorithm.h"
//Numerical solver for first-order explicit linear ODE, dx/dt+α(t)*x+f(t)=0
//Widely-accepted picking
//Variable step size
//DAEs? 
//Initial value problem? Boundary value problem (shooting)?
//Symplectic method

// GBDF (k=7,8) for Stiff
// Enright (k=1~4) for Stiff

vector<double> Scale(const vector<double>& v, double k) {
	vector<double> result(v.size());
	for (int i = 0; i < v.size(); ++i) {
		result[i] = v[i] * k;
	}
	return result;
}

vector<double> Add(const vector<double>& v1, const vector<double>& v2) {
    vector<double> result(v1.size());  
    for (int i = 0; i < v1.size(); ++i) {  
        result[i] = v1[i] + v2[i];  
    }  
    return result;
}

vector<string> Add(const vector<string>& v, string add) {
	vector<string> result(v.size());
	for (int i = 0; i < v.size(); ++i) {
		result[i] = v[i] + add;
	}
	return result;
}

vector<double> Substract(const vector<double>& v1, const vector<double>& v2) {
	vector<double> result(v1.size());
	for (int i = 0; i < v1.size(); ++i) {
		result[i] = v1[i] - v2[i];
	}
	return result;
}

vector<double> Square(const vector<double>& v) {
	vector<double> result(v.size());
	for (int i = 0; i < v.size(); ++i) {
		result[i] = v[i]*v[i];
	}
	return result;
}

vector<double> Multiply(const vector<double>& v1, const vector<double>& v2) {
	vector<double> result(v1.size());
	for (int i = 0; i < v1.size(); ++i) {
		result[i] = v1[i] * v2[i];
	}
	return result;
}

vector<double> Div(const vector<double>& v1, const vector<double>& v2) {
	vector<double> result(v1.size());
	for (int i = 0; i < v1.size(); ++i) {
		result[i] = v1[i] / v2[i];
	}
	return result;
}

double Max(const vector<double>& v) {
	double max = -999999;
	for (int i = 0; i < v.size(); ++i) {
		if (_Is_nan(abs(v[i])))
		{ 
			max = 1.0e20; 
		}
		else if (abs(v[i])>max)
		{
			max = abs(v[i]);
		}
	}

	return max;
}

double Tolerance(const vector<double>& v1, const vector<double>& v2) {
	double max=0,temp;
	for (int i = 0; i < v1.size(); ++i) {
		temp = abs(v1[i] - v2[i]);
		if (temp > max) {
			max = temp;
		}
	}
	return max;
}

//Euler_Forward V
vector<vector <double>> ODE_Solve_Euler_Forward(vector<double> (*Fun)(const vector<double> &x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i<len; ++i) {
		vector<double> F = Fun(x_solved[i - 1], time[i - 1], K);
		x_solved[i] = Add(Scale(F, h), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n+1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}
	
	return result;
}

//Euler Trapezoidal V(Improved Euler, with Predictor-Corrector Format)
vector<vector <double>> ODE_Solve_Euler_Trap(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> F = Fun(x_solved[i - 1], time[i - 1], K);
		x_solved[i] = Add(Scale(F, h), x_solved[i - 1]);
		F = Add(F, Fun(x_solved[i], time[i], K));
		x_solved[i] = Add(Scale(F, h / 2), x_solved[i - 1]);
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//Euler Backward V
vector<vector <double>> ODE_Solve_Euler_Backward(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
		vector<double> x_iter_n = Add(Scale(Fun(x_iter, time[i], K), h), x_solved[i - 1]);
		while (Tolerance(x_iter,x_iter_n)>1e-10) { // Fixed-Point Iteration
			x_iter = x_iter_n;
			x_iter_n = Add(Scale(Fun(x_iter, time[i], K), h), x_solved[i - 1]);
		}
		x_solved[i] = x_iter_n;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//Adamas-Bashforth Explicit V
vector<vector <double>> ODE_Solve_Adamas_Explicit(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) {  //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}
	
	if (k == 2) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 3 );
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -1);
			x_solved[i] = Add(Scale(Add(F1,F2), h/2), x_solved[i - 1]);
		}
	}
	else if (k == 3) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 23);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -16);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 5);
			x_solved[i] = Add(Scale(Add(Add(F1, F2),F3), h / 12), x_solved[i - 1]);
		}
	}
	else if (k == 4) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 55);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -59);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 37);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), - 9);
			x_solved[i] = Add(Scale(Add(Add(Add(F1,F2),F3),F4), h/24), x_solved[i - 1]);
		}
	}
	else if (k == 5) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 1901);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -2774);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 2616);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -1274);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 251);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(F1, F2), F3), F4), F5), h/720.0), x_solved[i - 1]);
		}
	}
	else if (k == 6) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 4277);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -7923);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 9982);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -7298);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 2877);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -475);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), h/1440.0), x_solved[i - 1]);
		}
	}
	else if (k == 7) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 198721.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -447288.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 705549.0);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -688256.0);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 407139.0);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -134472.0);
			vector<double> F7 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 19087.0);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7) , h / 60480.0), x_solved[i - 1]);
		}
	}
	else if (k == 8) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 434241.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -1152169.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 2183877.0);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -2664477.0);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 2102243.0);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -1041723.0);
			vector<double> F7 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 295767.0);
			vector<double> F8 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -36799.0);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), h / 120960.0), x_solved[i - 1]);
		}
	}
	else if (k == 9) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 14097247.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -43125206.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 95476786.0);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -139855262.0);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 137968480.0);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -91172642.0);
			vector<double> F7 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 38833486.0);
			vector<double> F8 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -9664106.0);
			vector<double> F9 = Scale(Fun(x_solved[i - 9], time[i - 9], K), 1070017.0);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), h / 3628800.0), x_solved[i - 1]);
		}
	}
	else if (k == 10) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 30277247.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -104995189.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 265932680.0);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -454661776.0);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 538363838.0);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -444772162.0);
			vector<double> F7 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 252618224.0);
			vector<double> F8 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -94307320.0);
			vector<double> F9 = Scale(Fun(x_solved[i - 9], time[i - 9], K), 20884811.0);
			vector<double> F10 = Scale(Fun(x_solved[i - 10], time[i - 10], K), -2082753.0);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10), h / 7257600.0), x_solved[i - 1]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//Adamas-Moulton Implicit V
vector<vector <double>> ODE_Solve_Adamas_Implicit(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) {  //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	if (k == 1) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 1), F1;
			
			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 1);
			vector<double> x_iter_n = Add(Scale(Add(F1, F2), h / 2), x_solved[i - 1]);
			
			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 1);
				x_iter_n = Add(Scale(Add(F1, F2), h / 2), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 2) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 8);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -1);
			
			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 5);
			vector<double> x_iter_n = Add(Scale(Add(Add(F1, F2), F3), h / 12), x_solved[i - 1]);
			
			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 5);
				x_iter_n = Add(Scale(Add(Add(F1, F2), F3), h / 12), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 3) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 19);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -5);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 1);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 9);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(F1, F2), F3), F4), h / 24), x_solved[i - 1]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 9);
				x_iter_n = Add(Scale(Add(Add(Add(F1, F2), F3), F4), h / 24), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 4) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 646);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -264);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 106);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -19);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 251);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(F1, F2), F3), F4), F5), h / 720), x_solved[i - 1]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 251);
				x_iter_n = Add(Scale(Add(Add(Add(Add(F1, F2), F3), F4), F5), h / 720), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 5) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 1427);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -798);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 482);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -173);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 27);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 475);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), h / 1440), x_solved[i - 1]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 475);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), h / 1440), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 6) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 65112);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -46461);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 37504);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -20211);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 6312);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -863);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 19087);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), h / 60480), x_solved[i - 1]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 19087);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), h / 60480), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 7) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 139849);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -121797);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 123133);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -88547);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 41499);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -11351);
			vector<double> F8 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 1375);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 36799);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), h / 120960), x_solved[i - 1]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 36799);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), h / 120960), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 8) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 4467094);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -4604594);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 5595358);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -5033120);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 3146338);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -1291214);
			vector<double> F8 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 312874);
			vector<double> F9 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -33953);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 1070017);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), h / 3628800), x_solved[i - 1]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 1070017);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), h / 3628800), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 9) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 9449717);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -11271304);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 16002320);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -17283646);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 13510082);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -7394032);
			vector<double> F8 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 2687864);
			vector<double> F9 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -583435);
			vector<double> F10 = Scale(Fun(x_solved[i - 9], time[i - 9], K), 57281);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 2082753);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10), h / 7257600), x_solved[i - 1]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 2082753);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10), h / 7257600), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 10) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 656185652);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -890175549);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 1446205080);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -1823311566);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 1710774528);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -1170597042);
			vector<double> F8 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 567450984);
			vector<double> F9 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -184776195);
			vector<double> F10 = Scale(Fun(x_solved[i - 9], time[i - 9], K), 36284876);
			vector<double> F11 = Scale(Fun(x_solved[i - 10], time[i - 10], K), -3250433);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 134211265);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10),F11), h / 479001600), x_solved[i - 1]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 134211265);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10), F11), h / 479001600), x_solved[i - 1]);
			}
			x_solved[i] = x_iter_n;
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//Adamas Hybrid PECE 4Exp-3Imp V
vector<vector <double>> ODE_Solve_Adamas_PECE(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int k = 4;
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) {  //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}
		
	for (int i = k; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> F1;
		vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 19);
		vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -5);
		vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 1);

		vector<double> F1_pre = Scale(Fun(x_solved[i - 1], time[i - 1], K), 55);
		vector<double> F2_pre = Scale(Fun(x_solved[i - 2], time[i - 2], K), -59);
		vector<double> F3_pre = Scale(Fun(x_solved[i - 3], time[i - 3], K), 37);
		vector<double> F4_pre = Scale(Fun(x_solved[i - 4], time[i - 4], K), -9);

		vector<double> x_pre = Add(Scale(Add(Add(Add(F1_pre, F2_pre), F3_pre), F4_pre), h / 24), x_solved[i - 1]);
		F1 = Scale(Fun(x_pre, time[i], K), 9);
		x_solved[i] = Add(Scale(Add(Add(Add(F1, F2), F3), F4), h / 24), x_solved[i - 1]);
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// Hamming Method  V
vector<vector <double>> ODE_Solve_Hamming(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int k = 4;
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) {  //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	for (int i = k; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> F1;
		vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 6);
		vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -3);
		
		vector<double> F1_pre = Scale(Fun(x_solved[i - 1], time[i - 1], K), 8);
		vector<double> F2_pre = Scale(Fun(x_solved[i - 2], time[i - 2], K), -4);
		vector<double> F3_pre = Scale(Fun(x_solved[i - 3], time[i - 3], K), 8);
		
		vector<double> x_pre = Add(Scale(Add(Add(F1_pre, F2_pre), F3_pre), h / 3.0), x_solved[i - 4]);
		F1 = Scale(Fun(x_pre, time[i], K), 3);
		x_solved[i] = Add(Scale(Add(Add(F1, F2), F3), h / 8.0), Add(Scale(x_solved[i - 1], 9.0 / 8.0), Scale(x_solved[i - 3], -1.0 / 8.0)));
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//Nystrom Explicit  V
vector<vector <double>> ODE_Solve_Nystrom_Explicit(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) {  //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	if (k == 2) {//Midpoint formula
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F = Fun(x_solved[i - 1], time[i - 1], K);
			x_solved[i] = Add(Scale(F, 2*h), x_solved[i - 2]);
		}
	}
	else if (k == 3) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 7);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -2);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 1);
			x_solved[i] = Add(Scale(Add(Add(F1, F2), F3), h / 3), x_solved[i - 2]);
		}
	}
	else if (k == 4) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 8);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -5);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 4);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -1);
			x_solved[i] = Add(Scale(Add(Add(Add(F1, F2), F3), F4), h / 3), x_solved[i - 2]);
		}
	}
	else if (k == 5) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 269);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -266);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 294);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -146);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 29);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(F1, F2), F3), F4), F5), h/90.0), x_solved[i - 2]);
		}
	}
	else if (k == 6) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 297.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -406.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 574.0);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -426.0);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 169.0);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -28.0);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), h / 90.0), x_solved[i - 2]);
		}
	}
	else if (k == 7) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 13613.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -23886.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 41193.0);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -40672.0);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 24183.0);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -8010.0);
			vector<double> F7 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 1139.0);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), h / 3780.0), x_solved[i - 2]);
		}
	}
	else if (k == 8) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 14720.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -31635.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 64440.0);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -79417.0);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 62928.0);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -31257.0);
			vector<double> F7 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 8888.0);
			vector<double> F8 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -1107.0);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), h / 3780.0), x_solved[i - 2]);
		}
	}
	else if (k == 9) { //doubt
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 439777.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -1208066.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 2839756.0);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -4195622.0);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 4155230.0);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -2750822.0);
			vector<double> F7 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 1173196.0);
			vector<double> F8 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -292226.0);
			vector<double> F9 = Scale(Fun(x_solved[i - 9], time[i - 9], K), 32377.0);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), h / 113400.0), x_solved[i - 2]);
		}
	}
	else if (k == 10) { //doubt
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 505625.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -1492898.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 3979084.0);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -6854054.0);
			vector<double> F5 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 8141878.0);
			vector<double> F6 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -6738470.0);
			vector<double> F7 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 3831628.0);
			vector<double> F8 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -1431554.0);
			vector<double> F9 = Scale(Fun(x_solved[i - 9], time[i - 9], K), 317209.0);
			vector<double> F10 = Scale(Fun(x_solved[i - 10], time[i - 10], K), -31648.0);
			x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10), h / 113400.0), x_solved[i - 2]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//Störmer Explicit
/*
vector<vector <double>> ODE_Solve_Stormer_Explicit(vector<double>(*Fun)(const vector<double>& x0, double t), double t_range[2], double step_size, vector<double> x0, int k) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) {  //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1]);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	if (k == 2) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F = Fun(x_solved[i - 1], time[i - 1]);
			x_solved[i] = Substract(Add(Scale(F, h * h), Scale(x_solved[i - 1], 2.0)), x_solved[i - 2]);
		}
	}
	else if (k == 3) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1]), 13.0);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2]), - 2.0);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3]), 1.0);
			x_solved[i] = Substract(Add(Scale(Add(Add(F1, F2), F3), h * h / 12.0), Scale(x_solved[i - 1], 2.0)), x_solved[i - 2]);
		}
	}
	else if (k == 4) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1 = Scale(Fun(x_solved[i - 1], time[i - 1]), 14);
			vector<double> F2 = Scale(Fun(x_solved[i - 2], time[i - 2]), -5);
			vector<double> F3 = Scale(Fun(x_solved[i - 3], time[i - 3]), 4);
			vector<double> F4 = Scale(Fun(x_solved[i - 4], time[i - 4]), -1);
			x_solved[i] = Substract(Add(Scale(Add(Add(Add(F1, F2), F3), F4), h * h / 12.0), Scale(x_solved[i - 1], 2.0)), x_solved[i - 2]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}
*/

//Milne-Simpson  V
vector<vector <double>> ODE_Solve_Milne_Simpson(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) { //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	if (k == 3) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 4);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), 1);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Scale(Add(Add(F1, F2), F3), h / 3), x_solved[i - 2]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Fun(x_iter, time[i], K);
				x_iter_n = Add(Scale(Add(Add(F1, F2), F3), h / 3), x_solved[i - 2]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 4) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 124);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), 24);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 4);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -1);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 29);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(F1, F2), F3), F4),F5), h / 90), x_solved[i - 2]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 29);
				x_iter_n = Add(Scale(Add(Add(Add(Add(F1, F2), F3), F4), F5), h / 90), x_solved[i - 2]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 5) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 129);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), 14);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 14);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -6);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 1);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 28);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), h / 90), x_solved[i - 2]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 28);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), h / 90), x_solved[i - 2]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 6) { //doubt
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 5640);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), 33);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 1328);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -87);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 264);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -37);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 1139);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), h / 3780), x_solved[i - 2]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 1139);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), h / 3780), x_solved[i - 2]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 7) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 5864);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -639);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 2448);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -1927);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 936);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -261);
			vector<double> F8 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 32);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 1107);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), h / 3780), x_solved[i - 2]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 1107);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), h / 3780), x_solved[i - 2]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 8) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 182584);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -42494);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 120088);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -116120);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 74728);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -31154);
			vector<double> F8 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 7624);
			vector<double> F9 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -833);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 32377);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), h / 113400), x_solved[i - 2]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 32377);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), h / 113400), x_solved[i - 2]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 9) {//doubt
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 189145);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), 68738);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 181324);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -207974);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 166582);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -92390);
			vector<double> F8 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 33868);
			vector<double> F9 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -7394);
			vector<double> F10 = Scale(Fun(x_solved[i - 9], time[i - 9], K), 729);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 31648);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10), h / 113400), x_solved[i - 2]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 31648);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10), h / 113400), x_solved[i - 2]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 10) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F1;
			vector<double> F2 = Scale(Fun(x_solved[i - 1], time[i - 1], K), 12908620);
			vector<double> F3 = Scale(Fun(x_solved[i - 2], time[i - 2], K), -6449433);
			vector<double> F4 = Scale(Fun(x_solved[i - 3], time[i - 3], K), 17067984);
			vector<double> F5 = Scale(Fun(x_solved[i - 4], time[i - 4], K), -22652334);
			vector<double> F6 = Scale(Fun(x_solved[i - 5], time[i - 5], K), 21705672);
			vector<double> F7 = Scale(Fun(x_solved[i - 6], time[i - 6], K), -15023790);
			vector<double> F8 = Scale(Fun(x_solved[i - 7], time[i - 7], K), 7335888);
			vector<double> F9 = Scale(Fun(x_solved[i - 8], time[i - 8], K), -2400729);
			vector<double> F10 = Scale(Fun(x_solved[i - 9], time[i - 9], K), 473164);
			vector<double> F11 = Scale(Fun(x_solved[i - 10], time[i - 10], K), -42505);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F1 = Scale(Fun(x_iter, time[i], K), 2046263);
			vector<double> x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10), F11), h / 7484400), x_solved[i - 2]);

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F1 = Scale(Fun(x_iter, time[i], K), 2046263);
				x_iter_n = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Add(Add(F1, F2), F3), F4), F5), F6), F7), F8), F9), F10), F11), h / 7484400), x_solved[i - 2]);
			}
			x_solved[i] = x_iter_n;
		}
	}
	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//R-K Middle point 2-order V
vector<vector <double>> ODE_Solve_RK_MP2(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		x_solved[i] = Add(Scale(K2, h), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//R-K Heun 2-order V
vector<vector <double>> ODE_Solve_RK_Heun2(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 2.0/3.0*h)), time[i - 1] + 2.0/3.0*h, K);
		x_solved[i] = Add(Scale(Add(K1,Scale(K2,3)), h / 4), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//R-K Heun 3-order V
vector<vector <double>> ODE_Solve_RK_Heun3(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h / 3)), time[i - 1] + h / 3, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 2.0 * h/3)), time[i - 1] + 2.0*h / 3, K);
		x_solved[i] = Add(Scale(Add(K1, Scale(K3, 3)), h / 4), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//R-K Kutta 3-order V
vector<vector <double>> ODE_Solve_RK_Kutta3(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h / 2)), time[i - 1] + h / 2, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, -h)),Scale(K2,2.0*h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(K1, Scale(K2, 4)),K3), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//R-K Nystrom 3-order V
vector<vector <double>> ODE_Solve_RK_Nystrom3(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 2.0/3.0*h)), time[i - 1] + 2.0/3.0*h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 2.0/3.0*h)), time[i - 1] + 2.0/3.0*h, K);
		x_solved[i] = Add(Scale(Add(Add(Scale(K1,2), Scale(K2, 3)), Scale(K3,3)), h / 8), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//R-K Kutta 4-order V
vector<vector <double>> ODE_Solve_RK_Kutta4(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h /3 )), time[i - 1] + h /3.0, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, - h / 3)), Scale(K2,h)), time[i - 1] + 2.0 * h / 3.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h)), Scale(K2, -h)), Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 3)), Scale(K3, 3)), K4), h / 8), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//R-K Classic 4-order V
vector<vector <double>> ODE_Solve_RK_C4(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1,Scale(K2,2)),Scale(K3,2)),K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//R-K Bogacki–Shampine  V
vector<vector <double>> ODE_Solve_RKBS(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.75 * h)), time[i - 1] + 0.75 * h, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, 2.0 / 9.0 * h)), Scale(K2, h / 3.0)), Scale(K3, h * 4.0 / 9.0)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(Scale(K1, 7), Scale(K2, 6)), Scale(K3, 8)), Scale(K4, 3)), h / 24), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//R-K Gill 4-order V
vector<vector <double>> ODE_Solve_RK_Gill4(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, (sqrt(2.0)-1.0)/2.0 * h)), Scale(K2,(1.0-sqrt(2.0)/2.0)*h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(Add(x_solved[i - 1], Scale(K2, -sqrt(2.0)/2.0*h)),Scale(K3, (1+sqrt(2.0)/2.0)*h)), time[i], K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2.0-sqrt(2.0))), Scale(K3, 2.0+sqrt(2.0))), K4), h / 6), x_solved[i - 1]);
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//RK 5 order Nystrom V
vector<vector <double>> ODE_Solve_RK_Nystrom5(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h /3.0)), time[i - 1] + h/3.0, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, 4.0/25.0*h)), Scale(K2, 6.0/25.0*h)), time[i - 1] + 2.0 * h / 5.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h/4.0)), Scale(K2, -3.0*h)), Scale(K3, 15.0/4.0*h)), time[i - 1] + h, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 6.0/81.0*h)), Scale(K2, 90.0/81.0*h)), Scale(K3, -50.0/81.0*h)), Scale(K4, 8.0/81.0*h)), time[i - 1] + 2.0/3.0*h, K);
		vector<double> K6 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 6.0 / 75.0*h)), Scale(K2, 36.0/75.0*h)), Scale(K3, 10.0/75.0*h)), Scale(K4, 8.0/75.0*h)), time[i - 1] + 4.0/5.0*h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(Scale(K1,23),Scale(K3,125)),Scale(K5,-81)),Scale(K6,125)), h / 192.0), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//RK 5 order Luther V
vector<vector <double>> ODE_Solve_RK_Luther5(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h)), time[i - 1] + h, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, h/2.0)), Scale(K2, h/2.0)), time[i - 1] + h, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, 14.0/64.0*h)), Scale(K2, 5.0/64.0*h)), Scale(K3, -3.0/64.0*h)), time[i - 1] + h/4.0, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -12.0 /96.0 *h)), Scale(K2, -12.0 /96.0 * h)), Scale(K3, 8.0/96.0*h)), Scale(K4, 64.0 / 96.0 * h)), time[i - 1] + h/2.0, K);
		vector<double> K6 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K2, -9.0 / 64.0 * h)), Scale(K3, 5.0/64.0*h)), Scale(K4, 16.0/64.0*h)), Scale(K5, 36.0/64.0*h)), time[i - 1] + 3.0/4.0*h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(Add(Scale(K1, 7), Scale(K3, 7)), Scale(K4, 32)), Scale(K5, 12)), Scale(K6,32)), h/90.0), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//RK 5 order Butcher V
vector<vector <double>> ODE_Solve_RK_Butcher5(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h/8.0)), time[i - 1] + h/8.0, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, h/4.0)), time[i - 1] + h/4.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h/2.0)), Scale(K2, -h)), Scale(K3, h)), time[i - 1] + h / 2.0, K);
		vector<double> K5 = Fun(Add(Add(x_solved[i - 1], Scale(K1, 3.0/16.0*h)), Scale(K4, 9.0 / 16.0 * h)), time[i - 1] + h *3.0/ 4.0, K);
		vector<double> K6 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -5.0/7.0 * h)), Scale(K2, 4.0/7.0*h)), Scale(K3, 12.0/7.0*h)), Scale(K4, -12.0/7.0*h)), Scale(K5, 8.0/7.0*h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(Add(Scale(K1, 7), Scale(K3, 32)), Scale(K4, 12)), Scale(K5, 32)), Scale(K6, 7)), h / 90.0), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//RK 5 order Butcher Method a  V
vector<vector <double>> ODE_Solve_RK_Butcher5a(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h / 5.0)), time[i - 1] + h / 5.0, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, h * 2.0 / 5.0)), time[i - 1] + h * 2.0 / 5.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 3.0 / 16.0)), Scale(K2, 0)), Scale(K3, h * 5.0 / 16.0)), time[i - 1] + h / 2.0, K);
		vector<double> K5 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, 1.0 / 4.0 * h)), Scale(K3, - 5.0 / 4.0 * h)), Scale(K4, 2.0 * h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Scale(K1, 1.0), Scale(K4, 4.0)), Scale(K5, 1.0)), h / 6.0), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//RK 5 order Butcher Method b  V
vector<vector <double>> ODE_Solve_RK_Butcher5b(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h / 5.0)), time[i - 1] + h / 5.0, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, h * 2.0 / 5.0)), time[i - 1] + h * 2.0 / 5.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 75.0 / 64.0)), Scale(K2, - 9.0/4.0 * h)), Scale(K3, h * 117.0 / 64.0)), time[i - 1] + h * 0.75, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -37.0 / 36.0 * h)), Scale(K2, 7.0 / 3.0 * h)), Scale(K3, -0.75 * h)), Scale(K4, 4.0 / 9.0 * h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(Scale(K1, 19.0 / 144.0), Scale(K3, 25.0 / 48.0)), Scale(K4, 2.0 / 9.0)), Scale(K5, 1.0 / 8.0)), h), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//RK 5 order Butcher Method c  V
vector<vector <double>> ODE_Solve_RK_Butcher5c(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h / 5.0)), time[i - 1] + h / 5.0, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, h * 2.0 / 5.0)), time[i - 1] + h * 2.0 / 5.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 161.0 / 192.0)), Scale(K2, -19.0 / 12.0 * h)), Scale(K3, h * 287.0 / 192.0)), time[i - 1] + h * 0.75, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -27.0 / 28.0 * h)), Scale(K2, 19.0 / 7.0 * h)), Scale(K3, -291.0 / 196.0 * h)), Scale(K4, 36.0 / 49.0 * h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(Scale(K1, 7.0 / 48.0), Scale(K3, 475.0 / 1008.0)), Scale(K4, 2.0 / 7.0)), Scale(K5, 7.0 / 72.0)), h), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K-F 45  V
vector<vector <double>> ODE_Solve_RKF45(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), h);
		vector<double> K2 = Scale(Fun(Add(x_solved[i - 1], Scale(K1, 1.0/4.0)), time[i - 1] + h / 4, K), h);
		vector<double> K3 = Scale(Fun(Add(Add(x_solved[i - 1], Scale(K1, 3.0/32.0)), Scale(K2, 9.0/32.0)), time[i - 1] + 3.0 * h / 8.0, K), h);
		vector<double> K4 = Scale(Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, 1932.0/2197.0)), Scale(K2, -7200.0/2197.0)), Scale(K3, 7296.0/2197.0)), time[i - 1] + 12.0*h/13.0, K), h);
		vector<double> K5 = Scale(Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 439.0 / 216.0)), Scale(K2, -8)), Scale(K3, 3680.0 / 513.0)), Scale(K4, -845.0/4104.0)), time[i - 1] + h, K), h);
		vector<double> K6 = Scale(Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, - 8.0 / 27.0)), Scale(K2, 2)), Scale(K3, -3544.0 / 2565.0)), Scale(K4, 1859.0 / 4104.0)), Scale(K5, -11.0/40.0)), time[i - 1] + h/2, K), h);
		x_solved[i] = Add(Add(Add(Add(Add(Scale(K1,16.0/135.0), Scale(K3, 6656.0/12825.0)), Scale(K4, 28561.0/56430.0)), Scale(K5,-9.0/50.0)), Scale(K6, 2.0/55.0)), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K-F 78  V
vector<vector <double>> ODE_Solve_RKF78(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Scale(Fun(x_solved[i - 1], time[i - 1], K), h);
		vector<double> K2 = Scale(Fun(Add(x_solved[i - 1], Scale(K1, 2.0 / 27.0)), time[i - 1] + h * 2.0 / 27.0, K), h);
		vector<double> K3 = Scale(Fun(Add(Add(x_solved[i - 1], Scale(K1, 1.0 / 36.0)), Scale(K2, 1.0 / 12.0)), time[i - 1] + h / 9.0, K), h);
		vector<double> K4 = Scale(Fun(Add(Add(x_solved[i - 1], Scale(K1, 1.0 / 24.0)), Scale(K3, 1.0 / 8.0)), time[i - 1] +  h / 6.0, K), h);
		vector<double> K5 = Scale(Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 5.0 / 12.0)), Scale(K2, 0)), Scale(K3, - 25.0 / 16.0)), Scale(K4, 25.0 / 16.0)), time[i - 1] + h*5.0/12.0, K), h);
		vector<double> K6 = Scale(Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 1.0 / 20.0)), Scale(K2, 0)), Scale(K3, 0)), Scale(K4, 1.0 / 4.0)), Scale(K5, 1.0 / 5.0)), time[i - 1] + h / 2, K), h);
		vector<double> K7 = Scale(Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -25.0 / 108.0)), Scale(K4, 125.0 / 108.0)), Scale(K5, - 65.0 / 27.0)), Scale(K6, 125.0 / 54.0)), time[i - 1] + h * 5.0 / 6.0, K), h);
		vector<double> K8 = Scale(Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 31.0 / 300.0)), Scale(K5, 61.0 / 225.0)), Scale(K6, -2.0 / 9.0)), Scale(K7, 13.0 / 900.0)), time[i - 1] + h / 6.0, K), h);
		vector<double> K9 = Scale(Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 2.0)), Scale(K4, -53.0 / 6.0)), Scale(K5, 704.0 / 45.0)), Scale(K6, -107.0 / 9.0)), Scale(K7, 67.0 / 90.0)), Scale(K8, 3.0)), time[i - 1] + h * 2.0 / 3.0, K), h);
		vector<double> K10 = Scale(Fun(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -91.0 / 108.0)), Scale(K4, 23.0 / 108.0)), Scale(K5, -976.0 / 135.0)), Scale(K6, 311.0 / 54.0)), Scale(K7, -19.0 / 60.0)), Scale(K8, 17.0 / 6.0)), Scale(K9, -1.0 / 12.0)), time[i - 1] + h * 1.0 / 3.0, K), h);
		vector<double> K11 = Scale(Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 2383.0 / 4100.0)), Scale(K4, -341.0 / 164.0)), Scale(K5, 4496.0 / 1025.0)), Scale(K6, -301.0 / 82.0)), Scale(K7, 2133.0 / 4100.0)), Scale(K8, 45.0 / 82.0)), Scale(K9, 45.0 / 164.0)), Scale(K10, 18.0 / 41.0)), time[i - 1] + h, K), h);
		vector<double> K12 = Scale(Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 3.0 / 205.0)), Scale(K6, -6.0 / 41.0)), Scale(K7, -3.0 / 205.0)), Scale(K8, -3.0 / 41.0)), Scale(K9, 3.0 / 41.0)), Scale(K10, 6.0 / 41.0)), time[i - 1], K), h);
		vector<double> K13 = Scale(Fun(Add(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -1777.0 / 4100.0)), Scale(K4, -341.0 / 164.0)), Scale(K5, 4496.0 / 1025.0)), Scale(K6, -289.0 / 82.0)), Scale(K7, 2193.0 / 4100.0)), Scale(K8, 51.0 / 82.0)), Scale(K9, 33.0 / 164.0)), Scale(K10, 12.0 / 41.0)), K12), time[i - 1] + h, K), h);
		x_solved[i] = Add(Add(Add(Add(Add(Add(Add(Scale(K6, 34.0 / 105.0), Scale(K7, 9.0 / 35.0)), Scale(K8, 9.0 / 35.0)), Scale(K9, 9.0 / 280.0)), Scale(K10, 9.0 / 280.0)), Scale(K12, 41.0/840.0)), Scale(K13, 41.0/840.0)), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K 6 order Butcher  V
vector<vector <double>> ODE_Solve_RK_Butcher6(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h/3.0)), time[i - 1] + h / 3.0, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 2.0 / 3.0*h)), time[i - 1] + 2.0 * h / 3.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h/12.0)), Scale(K2, h/3.0)), Scale(K3, -h/12.0)), time[i - 1] + h/3.0, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -h/16.0)), Scale(K2, 9.0/8.0*h)), Scale(K3, -3.0/16.0*h)), Scale(K4, -3.0/8.0*h)), time[i - 1] + 0.5*h, K);
		vector<double> K6 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K2, 9.0/8.0*h)), Scale(K3, -3.0/8.0*h)), Scale(K4, -3.0/4.0*h)), Scale(K5, 0.5*h)), time[i - 1] + h / 2.0, K);
		vector<double> K7 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 9.0/44.0 * h)), Scale(K2, -9.0 / 11.0 * h)), Scale(K3, 63.0/44.0 * h)), Scale(K4, 18.0/11.0 * h)), Scale(K6,-16.0/11.0*h)), time[i - 1] + h, K);
		x_solved[i] = Add(Add(Add(Add(Add(Add(Scale(K1, 11.0/120.0*h), Scale(K3, 27.0/40.0*h)), Scale(K4, 27.0/40.0*h)), Scale(K5, -4.0/15.0*h)), Scale(K6, -4.0/15.0*h)), Scale(K7,11.0/120.0*h)), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Dormand-Prince (ode45) V
vector<vector <double>> ODE_Solve_RKDP(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h / 5.0)), time[i - 1] + h / 5.0, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, 3.0 / 40.0 * h)), Scale(K2,9.0/40.0*h)), time[i - 1] + 3.0 * h / 10.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 44.0 / 45.0)), Scale(K2, -56.0* h / 15.0)), Scale(K3, 32.0* h / 9.0)), time[i - 1] + 4.0 * h / 5.0, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 19372.0 / 6561.0)), Scale(K2, -25360.0 / 2187.0 * h)), Scale(K3, 64448.0 / 6561.0 * h)), Scale(K4, - 212.0 / 729.0 * h)), time[i - 1] + 8.0 / 9.0 * h, K);
		vector<double> K6 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 9017.0 / 3168.0 * h)), Scale(K2, -355.0 / 33.0 * h)), Scale(K3, 46732.0 / 5247.0 * h)), Scale(K4, 49.0/176.0 * h)), Scale(K5, - 5103.0 / 18656.0 * h)), time[i - 1] + h, K);
		vector<double> K7 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 35.0 / 384.0 * h)), Scale(K3, 500.0 / 1113.0 * h)), Scale(K4, 125.0 / 192.0 * h)), Scale(K5, - 2187.0 / 6784.0 * h)), Scale(K6, 11.0 / 84.0 * h)), time[i - 1] + h, K);
		x_solved[i] = Add(Add(Add(Add(Add(Add(Scale(K1, 5179.0 / 57600.0 * h), Scale(K3, 7571.0 / 16695.0 * h)), Scale(K4, 393.0 / 640.0 * h)), Scale(K5, -92097.0 / 339200.0 * h)), Scale(K6, 187.0 / 2100.0 * h)),Scale(K7, h/40.0) ), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Cash-Karp  V 
vector<vector <double>> ODE_Solve_RKCK(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h / 5.0)), time[i - 1] + h / 5.0, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, 3.0 / 40.0 * h)), Scale(K2, 9.0 / 40.0 * h)), time[i - 1] + 3.0 * h / 10.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 3.0 / 10.0)), Scale(K2, -9.0 * h / 10.0)), Scale(K3, 6.0 * h / 5.0)), time[i - 1] + 3.0 * h / 5.0, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, h * -11.0/54.0)), Scale(K2, 5.0 / 2.0 * h)), Scale(K3, - 70.0 / 27.0 * h)), Scale(K4, 35.0 / 27.0 * h)), time[i - 1] +  h, K);
		vector<double> K6 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 1631.0 / 55296.0 * h)), Scale(K2, 175.0 / 512.0 * h)), Scale(K3, 575.0 / 13828.0 * h)), Scale(K4, 44275.0 / 110592.0 * h)), Scale(K5, 253.0 / 4096.0 * h)), time[i - 1] + 7.0 / 8.0 * h, K);
		x_solved[i] = Add(Add(Add(Add(Add(Scale(K1, 2825.0 / 27648.0 * h), Scale(K3, 18575.0 / 48384.0 * h)), Scale(K4, 13525.0 / 55296.0 * h)), Scale(K5, 277.0 / 14336.0 * h)), Scale(K6, 1.0 / 4.0 * h)), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K HIHA  V
vector<vector <double>> ODE_Solve_RKHH(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h * 2.0 / 9.0)), time[i - 1] + h * 2.0 / 9.0, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, h / 12.0)), Scale(K2, h / 4.0)), time[i - 1] + h / 3.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 1.0 / 8.0)), Scale(K2, 0)), Scale(K3, 3.0 * h / 8.0)), time[i - 1] + 0.5 * h, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 91.0 / 500.0)), Scale(K2, -27.0 / 100.0 * h)), Scale(K3, 78.0 / 125.0 * h)), Scale(K4, 8.0 / 125.0 * h)), time[i - 1] + 3.0 / 5.0 * h, K);
		vector<double> K6 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -11.0 / 20.0 * h)), Scale(K2, 27.0 / 20.0 * h)), Scale(K3, h * 12.0 / 5.0)), Scale(K4, - 36.0 / 5.0 * h)), Scale(K5, 5.0 * h)), time[i - 1] + h, K);
		vector<double> K7 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 1.0 / 12.0 * h)), Scale(K3, 27.0 / 32.0 * h)), Scale(K4, -4.0 / 3.0 * h)), Scale(K5, 125.0 / 96.0 * h)), Scale(K6, 5.0 / 48.0 * h)), time[i - 1] + h, K);
		x_solved[i] = Add(Add(Add(Add(Add(Add(Scale(K1, 2.0 / 15.0 * h), Scale(K3, 27.0 / 80.0 * h)), Scale(K4, -2.0 / 15.0 * h)), Scale(K5, 25.0 / 48.0 * h)), Scale(K6, 1.0 / 24.0 * h)), Scale(K7, h / 10.0)), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K 54 7S V
vector<vector <double>> ODE_Solve_RK7S(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h * 2.0 / 9.0)), time[i - 1] + h * 2.0 / 9.0, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, h /12.0)), Scale(K2, h / 4.0)), time[i - 1] + h / 3.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 55.0 / 324.0)), Scale(K2, - 25.0 * h / 108.0)), Scale(K3, 50.0 * h / 81.0)), time[i - 1] + 5.0 * h / 9.0, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 83.0/330.0)), Scale(K2, -13.0/22.0 * h)), Scale(K3, 61.0/66.0 * h)), Scale(K4, 9.0 / 110.0 * h)), time[i - 1] + 2.0 / 3.0 * h, K);
		vector<double> K6 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, - 19.0 / 28.0 * h)), Scale(K2, 9.0 / 4.0 * h)), Scale(K3, h / 7.0)), Scale(K4, - 27.0 / 7.0 * h)), Scale(K5, 22.0 / 7.0 * h)), time[i - 1] + h, K);
		vector<double> K7 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, 19.0 / 200.0 * h)), Scale(K3, 3.0 / 5.0 * h)), Scale(K4, - 243.0 / 400.0 * h)), Scale(K5, 33.0 / 40.0 * h)), Scale(K6, 7.0 / 80.0 * h)), time[i - 1] + h, K);
		x_solved[i] = Add(Add(Add(Add(Add(Add(Scale(K1, 431.0 / 5000.0 * h), Scale(K3, 333.0 / 500.0 * h)), Scale(K4, - 7857.0 / 10000.0 * h)), Scale(K5, 957.0 / 1000.0 * h)), Scale(K6, 193.0 / 2000.0 * h)), Scale(K7, - h / 50.0)), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K 54 6M V
vector<vector <double>> ODE_Solve_RK6M(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, h / 5.0)), time[i - 1] + h / 5.0, K);
		vector<double> K3 = Fun(Add(Add(x_solved[i - 1], Scale(K1, h * 3.0 / 40.0)), Scale(K2, h * 9.0 / 40.0)), time[i - 1] + h * 3.0 / 10.0, K);
		vector<double> K4 = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 3.0 / 10.0)), Scale(K2, - 9.0 * h / 10.0)), Scale(K3, 6.0 * h / 5.0)), time[i - 1] + 3.0 * h / 5.0, K);
		vector<double> K5 = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, h * 226.0 / 729.0)), Scale(K2, -25.0 / 27.0 * h)), Scale(K3, 880.0 / 729.0 * h)), Scale(K4, 55.0 / 729.0 * h)), time[i - 1] + 2.0 / 3.0 * h, K);
		vector<double> K6 = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1, -181.0 / 270.0 * h)), Scale(K2, 5.0 / 2.0 * h)), Scale(K3, - 266.0 * h / 297.0)), Scale(K4, - 91.0 / 27.0 * h)), Scale(K5, 189.0 / 55.0 * h)), time[i - 1] + h, K);
		x_solved[i] = Add(Add(Add(Add(Add(Scale(K1, 31.0 / 540.0 * h), Scale(K3, 190.0 / 297.0 * h)), Scale(K4, -145.0 / 108.0 * h)), Scale(K5, 351.0 / 220.0 * h)), Scale(K6, h / 20.0)), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//BDF k=2~6 adaptive for Stiff prob V
vector<vector <double>> ODE_Solve_BDF(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) { //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	if (k == 2) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			
			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Scale(x_solved[i-1],4.0/3.0), Scale(x_solved[i-2],-1.0/3.0)),Scale(F,2.0/3.0*h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Scale(x_solved[i - 1], 4.0 / 3.0), Scale(x_solved[i - 2], -1.0 / 3.0)), Scale(F, 2.0 / 3.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 3) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Scale(x_solved[i - 1], 3.0* 6.0 / 11.0), Scale(x_solved[i - 2], -3.0 / 2.0*6.0/11.0)), Scale(x_solved[i-3],2.0/11.0)), Scale(F, 6.0 / 11.0 * h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Scale(x_solved[i - 1], 3.0 * 6.0 / 11.0), Scale(x_solved[i - 2], -3.0 / 2.0 * 6.0 / 11.0)), Scale(x_solved[i - 3], 2.0 / 11.0)), Scale(F, 6.0 / 11.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 4) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Add(Scale(x_solved[i - 1], 4.0 * 12.0 / 25.0), Scale(x_solved[i - 2], -3.0* 12.0 / 25.0)), Scale(x_solved[i - 3], 4.0/3.0*12.0/25.0)), Scale(x_solved[i-4],-1.0/4.0*12.0/25.0)), Scale(F, 12.0/25.0 * h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Add(Scale(x_solved[i - 1], 4.0 * 12.0 / 25.0), Scale(x_solved[i - 2], -3.0 * 12.0 / 25.0)), Scale(x_solved[i - 3], 4.0 / 3.0 * 12.0 / 25.0)), Scale(x_solved[i - 4], -1.0 / 4.0 * 12.0 / 25.0)), Scale(F, 12.0 / 25.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 5) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Add(Add(Scale(x_solved[i - 1], 5.0 * 60.0/137.0), Scale(x_solved[i - 2], -5.0 * 60.0/137.0)), Scale(x_solved[i - 3], 10.0 / 3.0 * 60.0/137.0)), Scale(x_solved[i - 4], -5.0 / 4.0 * 60.0/137.0)), Scale(x_solved[i-5],1.0/5.0*60.0/137.0)), Scale(F, 60.0 / 137.0 * h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Add(Add(Scale(x_solved[i - 1], 5.0 * 60.0 / 137.0), Scale(x_solved[i - 2], -5.0 * 60.0 / 137.0)), Scale(x_solved[i - 3], 10.0 / 3.0 * 60.0 / 137.0)), Scale(x_solved[i - 4], -5.0 / 4.0 * 60.0 / 137.0)), Scale(x_solved[i - 5], 1.0 / 5.0 * 60.0 / 137.0)), Scale(F, 60.0 / 137.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 6) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Add(Add(Add(Scale(x_solved[i - 1], 6.0 * 60.0 / 147.0), Scale(x_solved[i - 2], -15.0/2.0 * 60.0 / 147.0)), Scale(x_solved[i - 3], 20.0 / 3.0 * 60.0 / 147.0)), Scale(x_solved[i - 4], -15.0 / 4.0 * 60.0 / 147.0)), Scale(x_solved[i - 5], 6.0 / 5.0 * 60.0 / 147.0)), Scale(x_solved[i-6],-1.0/6.0*60.0/147.0)), Scale(F, 60.0 / 147.0 * h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Add(Add(Add(Scale(x_solved[i - 1], 6.0 * 60.0 / 147.0), Scale(x_solved[i - 2], -15.0 / 2.0 * 60.0 / 147.0)), Scale(x_solved[i - 3], 20.0 / 3.0 * 60.0 / 147.0)), Scale(x_solved[i - 4], -15.0 / 4.0 * 60.0 / 147.0)), Scale(x_solved[i - 5], 6.0 / 5.0 * 60.0 / 147.0)), Scale(x_solved[i - 6], -1.0 / 6.0 * 60.0 / 147.0)), Scale(F, 60.0 / 147.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//NDF  V
vector<vector <double>> ODE_Solve_NDF(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) { //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	if (k == 2) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(x_solved[i - 1], Scale(F, h)), Scale(Add(Add(x_iter, Scale(x_solved[i - 1], -2.0)), x_solved[i - 2]), - 37.0 / 200.0));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(x_solved[i - 1], Scale(F, h)), Scale(Add(Add(x_iter, Scale(x_solved[i - 1], -2.0)), x_solved[i - 2]), -37.0 / 200.0));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 3) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Scale(x_solved[i - 1], 4.0 / 3.0), Scale(x_solved[i - 2], -1.0 / 3.0)), Scale(F, 2.0 / 3.0 * h)), Scale(Substract(Add(Add(x_iter, Scale(x_solved[i - 1], -3.0)), Scale(x_solved[i - 2], 3.0)), x_solved[i - 3]), -1.0 / 9.0));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Scale(x_solved[i - 1], 4.0 / 3.0), Scale(x_solved[i - 2], -1.0 / 3.0)), Scale(F, 2.0 / 3.0 * h)), Scale(Substract(Add(Add(x_iter, Scale(x_solved[i - 1], -3.0)), Scale(x_solved[i - 2], 3.0)), x_solved[i - 3]), -1.0 / 9.0));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 4) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Add(Scale(x_solved[i - 1], 18.0 / 11.0), Scale(x_solved[i - 2], -9.0 / 11.0)), Scale(x_solved[i - 3], 2.0 / 11.0)), Scale(F, 6.0 / 11.0 * h)), Scale(Add(Add(Add(Add(x_iter, Scale(x_solved[i - 1], -4.0)), Scale(x_solved[i - 2], 6.0)), Scale(x_solved[i - 3], -4.0)), x_solved[i - 4]), 6.0 / 11.0 * -0.0823 * (1.0 + 1.0 / 2.0 + 1.0 / 3.0)));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Add(Scale(x_solved[i - 1], 18.0 / 11.0), Scale(x_solved[i - 2], -9.0 / 11.0)), Scale(x_solved[i - 3], 2.0 / 11.0)), Scale(F, 6.0 / 11.0 * h)), Scale(Add(Add(Add(Add(x_iter, Scale(x_solved[i - 1], -4.0)), Scale(x_solved[i - 2], 6.0)), Scale(x_solved[i - 3], -4.0)), x_solved[i - 4]), 6.0 / 11.0 * -0.0823 * (1.0 + 1.0 / 2.0 + 1.0 / 3.0)));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 5) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Add(Add(Scale(x_solved[i - 1], 48.0 / 25.0), Scale(x_solved[i - 2], -36.0 / 25.0)), Scale(x_solved[i - 3], 16.0 / 25.0)), Scale(x_solved[i - 4], -3.0 / 25.0)), Scale(F, 12.0 / 25.0 * h)), Scale(Substract(Add(Add(Add(Add(x_iter, Scale(x_solved[i - 1], -5.0)), Scale(x_solved[i - 2], 10.0)), Scale(x_solved[i - 3], -10.0)), Scale(x_solved[i - 4], 5.0)), x_solved[i - 5]), 12.0 / 25.0 * -0.0415 * (1.0 + 1.0 / 2.0 + 1.0 / 3.0 + 1.0 / 4.0)));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Add(Add(Scale(x_solved[i - 1], 48.0 / 25.0), Scale(x_solved[i - 2], -36.0 / 25.0)), Scale(x_solved[i - 3], 16.0 / 25.0)), Scale(x_solved[i - 4], -3.0 / 25.0)), Scale(F, 12.0 / 25.0 * h)), Scale(Substract(Add(Add(Add(Add(x_iter, Scale(x_solved[i - 1], -5.0)), Scale(x_solved[i - 2], 10.0)), Scale(x_solved[i - 3], -10.0)), Scale(x_solved[i - 4], 5.0)), x_solved[i - 5]), 12.0 / 25.0 * -0.0415 * (1.0 + 1.0 / 2.0 + 1.0 / 3.0 + 1.0 / 4.0)));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 6) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			vector<double> F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Add(Add(Scale(x_solved[i - 1], 300.0 / 137.0), Scale(x_solved[i - 2], -300.0 / 137.0)), Scale(x_solved[i - 3], 200.0 / 137.0)), Scale(x_solved[i - 4], -75.0 / 137.0)), Scale(x_solved[i - 5], 12.0 / 137.0)), Scale(F, 60.0 / 137.0 * h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Add(Add(Scale(x_solved[i - 1], 300.0 / 137.0), Scale(x_solved[i - 2], -300.0 / 137.0)), Scale(x_solved[i - 3], 200.0 / 137.0)), Scale(x_solved[i - 4], -75.0 / 137.0)), Scale(x_solved[i - 5], 12.0 / 137.0)), Scale(F, 60.0 / 137.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// Runge-Kutta Implicit (Gauss-Legendre, s=1~3) for Stiff V
vector<vector <double>> ODE_Solve_RK_Implicit(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	if (s == 1) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1] + h / 2, K);
			vector<double> K1_iter_n=Fun(Add(x_solved[i-1],Scale(K1_iter,h/2)),time[i-1]+h/2, K);
			while (Tolerance(K1_iter, K1_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K1_iter_n = Fun(Add(x_solved[i - 1], Scale(K1_iter, h / 2)), time[i - 1] + h / 2, K);
			}
			
			x_solved[i] = Add(Scale(K1_iter_n, h), x_solved[i - 1]);
		}
	}
	else if (s == 2) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1] + (1.0 / 2.0 - sqrt(3.0) / 6.0) * h, K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (1.0 / 2.0 + sqrt(3.0) / 6.0) * h, K);

			vector<double> K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 4)), Scale(K2_iter, (1.0 / 4.0 - sqrt(3.0) / 6.0) * h)), time[i - 1] + (1.0 / 2.0 - sqrt(3.0) / 6.0) * h, K);
			vector<double> K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K2_iter, h / 4)), Scale(K1_iter, (1.0 / 4.0 + sqrt(3.0) / 6.0) * h)), time[i - 1] + (1.0 / 2.0 + sqrt(3.0) / 6.0) * h, K);
			
			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 4)), Scale(K2_iter, (1.0 / 4.0 - sqrt(3.0) / 6.0) * h)), time[i - 1] + (1.0 / 2.0 - sqrt(3.0) / 6.0) * h, K);
				K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K2_iter, h / 4)), Scale(K1_iter, (1.0 / 4.0 + sqrt(3.0) / 6.0) * h)), time[i - 1] + (1.0 / 2.0 + sqrt(3.0) / 6.0) * h, K);
			}

			x_solved[i] = Add(Scale(Add(K1_iter_n,K2_iter_n), h/2), x_solved[i - 1]);
		}
	}
	else if (s == 3) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1] + (1.0 / 2.0 - sqrt(15.0) / 10.0) * h, K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 1.0 / 2.0 * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + (1.0 / 2.0 + sqrt(15.0) / 10.0) * h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h *5.0/36.0)), Scale(K2_iter, (2.0 / 9.0 - sqrt(15.0) / 15.0) * h)), Scale(K3_iter, (5.0/36.0-sqrt(15.0)/30.0)*h)), time[i - 1] + (1.0 / 2.0 - sqrt(15.0) / 10.0) * h, K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (5.0 / 36.0+sqrt(15.0)/24.0))), Scale(K2_iter, 2.0 / 9.0 * h)), Scale(K3_iter, (5.0 / 36.0 - sqrt(15.0) / 24.0) * h)), time[i - 1] + 1.0 / 2.0 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (5.0 / 36.0+sqrt(15.0)/30.0))), Scale(K2_iter, (2.0 / 9.0 + sqrt(15.0) / 15.0) * h)), Scale(K3_iter, 5.0 / 36.0 * h)), time[i - 1] + (1.0 / 2.0 + sqrt(15.0) / 10.0) * h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 5.0 / 36.0)), Scale(K2_iter, (2.0 / 9.0 - sqrt(15.0) / 15.0) * h)), Scale(K3_iter, (5.0 / 36.0 - sqrt(15.0) / 30.0) * h)), time[i - 1] + (1.0 / 2.0 - sqrt(15.0) / 10.0) * h, K);
				K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (5.0 / 36.0 + sqrt(15.0) / 24.0))), Scale(K2_iter, 2.0 / 9.0 * h)), Scale(K3_iter, (5.0 / 36.0 - sqrt(15.0) / 24.0) * h)), time[i - 1] + 1.0 / 2.0 * h, K);
				K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (5.0 / 36.0 + sqrt(15.0) / 30.0))), Scale(K2_iter, (2.0 / 9.0 + sqrt(15.0) / 15.0) * h)), Scale(K3_iter, 5.0 / 36.0 * h)), time[i - 1] + (1.0 / 2.0 + sqrt(15.0) / 10.0) * h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Scale(K1_iter_n,5), Scale(K2_iter_n,8)),Scale(K3_iter_n,5)), h / 18), x_solved[i - 1]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

//Gear K=2~6 V
vector<vector <double>> ODE_Solve_Gear(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int k, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < k; ++i) { //R-K Classic 4-order V
		vector<double> K1 = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2 = Fun(Add(x_solved[i - 1], Scale(K1, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K3 = Fun(Add(x_solved[i - 1], Scale(K2, 0.5 * h)), time[i - 1] + 0.5 * h, K);
		vector<double> K4 = Fun(Add(x_solved[i - 1], Scale(K3, h)), time[i - 1] + h, K);
		x_solved[i] = Add(Scale(Add(Add(Add(K1, Scale(K2, 2)), Scale(K3, 2)), K4), h / 6), x_solved[i - 1]);
		time[i] = time[i - 1] + h;
	}

	if (k == 2) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F;
			vector<double> X1 = Scale(x_solved[i-2], -1.0/3.0);
			vector<double> X2 = Scale(x_solved[i-1], 4.0/3.0);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(X1,X2),Scale(F,2.0/3.0*h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(X1, X2), Scale(F, 2.0 / 3.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 3) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F;
			vector<double> X1 = Scale(x_solved[i - 3], 2.0/11.0);
			vector<double> X2 = Scale(x_solved[i - 2], -9.0/11.0);
			vector<double> X3 = Scale(x_solved[i - 1], 18.0/11.0);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(X1, X2), X3), Scale(F, 6.0/11.0*h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(X1, X2), X3), Scale(F, 6.0 / 11.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 4) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F;
			vector<double> X1 = Scale(x_solved[i - 4], -3.0/25.0);
			vector<double> X2 = Scale(x_solved[i - 3], 16.0/25.0);
			vector<double> X3 = Scale(x_solved[i - 2], -36.0/25.0);
			vector<double> X4 = Scale(x_solved[i - 1], 48.0/25.0);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Add(X1, X2), X3), X4), Scale(F, 12.0/25.0 * h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Add(X1, X2), X3), X4), Scale(F, 12.0 / 25.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 5) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F;
			vector<double> X1 = Scale(x_solved[i - 5], 12.0 / 137.0);
			vector<double> X2 = Scale(x_solved[i - 4], -75.0 / 137.0);
			vector<double> X3 = Scale(x_solved[i - 3], 200.0/137.0);
			vector<double> X4 = Scale(x_solved[i - 2], -300.0/137.0);
			vector<double> X5 = Scale(x_solved[i - 1], 300.0 / 137.0);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Add(Add(X1, X2), X3), X4), X5), Scale(F, 60.0/137.0 * h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Add(Add(X1, X2), X3), X4), X5), Scale(F, 60.0 / 137.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}
	else if (k == 6) {
		for (int i = k; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> F;
			vector<double> X1 = Scale(x_solved[i - 6], -10.0/147.0);
			vector<double> X2 = Scale(x_solved[i - 5], 72.0/147.0);
			vector<double> X3 = Scale(x_solved[i - 4], -225.0/147.0);
			vector<double> X4 = Scale(x_solved[i - 3], 400.0/147.0);
			vector<double> X5 = Scale(x_solved[i - 2], -450.0/147.0);
			vector<double> X6 = Scale(x_solved[i - 1], 360.0/147.0);

			vector<double> x_iter = Add(Scale(Fun(x_solved[i - 1], time[i - 1], K), h), x_solved[i - 1]);
			F = Fun(x_iter, time[i], K);
			vector<double> x_iter_n = Add(Add(Add(Add(Add(Add(X1, X2), X3), X4), X5), X6), Scale(F, 60.0/147.0 * h));

			while (Tolerance(x_iter, x_iter_n) > 1e-10) { // Fixed-Point Iteration
				x_iter = x_iter_n;
				F = Fun(x_iter, time[i], K);
				x_iter_n = Add(Add(Add(Add(Add(Add(X1, X2), X3), X4), X5), X6), Scale(F, 60.0 / 147.0 * h));
			}
			x_solved[i] = x_iter_n;
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Implicit RadauIA  V
vector<vector <double>> ODE_Solve_RK_RIA(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	if (s == 2) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + h * 2.0 / 3.0, K);

			vector<double> K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 4.0)), Scale(K2_iter, - h / 4.0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 4.0)), Scale(K2_iter, h * 5.0 / 12.0)), time[i - 1] + h * 2.0 / 3.0, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 4.0)), Scale(K2_iter, -h / 4.0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 4.0)), Scale(K2_iter, h * 5.0 / 12.0)), time[i - 1] + h * 2.0 / 3.0, K);
			}

			x_solved[i] = Add(Add(Scale(K1_iter_n, 1.0 / 4.0 * h), Scale(K2_iter_n, 3.0 / 4.0 * h)), x_solved[i - 1]);
		}
	}
	else if (s == 3) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (6.0 - sqrt(6.0)) / 10.0 * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + (6.0 + sqrt(6.0)) / 10.0 * h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 9.0)), Scale(K2_iter, (-1.0 - sqrt(6.0)) / 18.0 * h)), Scale(K3_iter, (-1.0 + sqrt(6.0)) / 18.0 * h)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 9.0)), Scale(K2_iter, (88.0 + 7.0 * sqrt(6.0)) / 360.0 * h)), Scale(K3_iter, (88.0 - 43.0 * sqrt(6.0)) / 360.0 * h)), time[i - 1] + (6.0 - sqrt(6.0)) / 10.0 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 9.0)), Scale(K2_iter, (88.0 + 43.0 * sqrt(6.0)) / 360.0 * h)), Scale(K3_iter,(88.0 - 7.0 * sqrt(6.0)) / 360.0 * h)), time[i - 1] + (6.0 + sqrt(6.0)) / 10.0 * h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 9.0)), Scale(K2_iter, (-1.0 - sqrt(6.0)) / 18.0 * h)), Scale(K3_iter, (-1.0 + sqrt(6.0)) / 18.0 * h)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 9.0)), Scale(K2_iter, (88.0 + 7.0 * sqrt(6.0)) / 360.0 * h)), Scale(K3_iter, (88.0 - 43.0 * sqrt(6.0)) / 360.0 * h)), time[i - 1] + (6.0 - sqrt(6.0)) / 10.0 * h, K);
				K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 9.0)), Scale(K2_iter, (88.0 + 43.0 * sqrt(6.0)) / 360.0 * h)), Scale(K3_iter, (88.0 - 7.0 * sqrt(6.0)) / 360.0 * h)), time[i - 1] + (6.0 + sqrt(6.0)) / 10.0 * h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Scale(K1_iter_n, 1.0 / 9.0), Scale(K2_iter_n, (16.0 + sqrt(6.0)) / 36.0)), Scale(K3_iter_n, (16.0 - sqrt(6.0)) / 36.0)), h), x_solved[i - 1]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Implicit RadauIIA  V
vector<vector <double>> ODE_Solve_RK_RIIA(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	if (s == 2) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1] + h / 3.0, K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 5.0 / 12.0)), Scale(K2_iter, -h / 12.0)), time[i - 1] + h / 3.0, K);
			vector<double> K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 3.0 / 4.0)), Scale(K2_iter, h / 4.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 5.0 / 12.0)), Scale(K2_iter, -h / 12.0)), time[i - 1] + h / 3.0, K);
				K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 3.0 / 4.0)), Scale(K2_iter, h / 4.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Add(Scale(K1_iter_n, 3.0 / 4.0 * h), Scale(K2_iter_n, 1.0 / 4.0 * h)), x_solved[i - 1]);
		}
	}
	else if (s == 3) {//Stiff
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1] + (4.0 - sqrt(6.0)) / 10.0 * h, K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (4.0 + sqrt(6.0)) / 10.0 * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (88.0 - 7.0 * sqrt(6.0)) / 360.0)), Scale(K2_iter, (296.0 - 169.0 * sqrt(6.0)) / 1800.0 * h)), Scale(K3_iter, (-2.0 + 3.0 * sqrt(6.0)) / 225.0 * h)), time[i - 1] + (4.0 - sqrt(6.0)) / 10.0 * h, K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (296.0 + 169.0 * sqrt(6.0)) / 1800.0)), Scale(K2_iter, (88.0 + 7.0 * sqrt(6.0)) / 360.0 * h)), Scale(K3_iter, (-2.0 - 3.0 * sqrt(6.0)) / 225.0 * h)), time[i - 1] + (4.0 + sqrt(6.0)) / 10.0 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (16.0 - sqrt(6.0)) / 36.0)), Scale(K2_iter, (16.0 + sqrt(6.0)) / 36.0 * h)), Scale(K3_iter, h / 9.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (88.0 - 7.0 * sqrt(6.0)) / 360.0)), Scale(K2_iter, (296.0 - 169.0 * sqrt(6.0)) / 1800.0 * h)), Scale(K3_iter, (-2.0 + 3.0 * sqrt(6.0)) / 225.0 * h)), time[i - 1] + (4.0 - sqrt(6.0)) / 10.0 * h, K);
				K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (296.0 + 169.0 * sqrt(6.0)) / 1800.0)), Scale(K2_iter, (88.0 + 7.0 * sqrt(6.0)) / 360.0 * h)), Scale(K3_iter, (-2.0 - 3.0 * sqrt(6.0)) / 225.0 * h)), time[i - 1] + (4.0 + sqrt(6.0)) / 10.0 * h, K);
				K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (16.0 - sqrt(6.0)) / 36.0)), Scale(K2_iter, (16.0 + sqrt(6.0)) / 36.0 * h)), Scale(K3_iter, h / 9.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Scale(K1_iter_n, (16.0 - sqrt(6.0)) / 36.0), Scale(K2_iter_n, (16.0 + sqrt(6.0)) / 36.0)), Scale(K3_iter_n, 1.0 / 9.0)), h), x_solved[i - 1]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Lobatto IIIA  V
vector<vector <double>> ODE_Solve_RK_LIIIA(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	if (s == 2) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 2.0)), Scale(K2_iter, h / 2.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 2.0)), Scale(K2_iter, h / 2.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Add(Scale(K1_iter_n, h / 2.0), Scale(K2_iter_n, h / 2.0)), x_solved[i - 1]);
		}
	}
	else if (s == 3) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 5.0 / 24.0)), Scale(K2_iter, h / 3.0)), Scale(K3_iter, - h / 24.0)), time[i - 1] + 0.5 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 2.0 * h / 3.0)), Scale(K3_iter, h / 6.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 5.0 / 24.0)), Scale(K2_iter, h / 3.0)), Scale(K3_iter, -h / 24.0)), time[i - 1] + 0.5 * h, K);
				K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 2.0 * h / 3.0)), Scale(K3_iter, h / 6.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Scale(K1_iter_n, 1.0 / 6.0), Scale(K2_iter_n, 2.0 / 3.0)), Scale(K3_iter_n, 1.0 / 6.0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 4) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, (11.0 + sqrt(5.0)) / 120.0 * h)), Scale(K2_iter, (25.0 - sqrt(5.0)) / 120.0 * h)), Scale(K3_iter, (25.0 - 13.0 * sqrt(5.0)) / 120.0 * h)), Scale(K4_iter, (-1.0 + sqrt(5.0)) / 120.0 * h)), time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, (11.0 - sqrt(5.0)) / 120.0 * h)), Scale(K2_iter, (25.0 + 13.0 * sqrt(5.0)) / 120.0 * h)), Scale(K3_iter, (25.0 + sqrt(5.0)) / 120.0 * h)), Scale(K4_iter, (-1.0 - sqrt(5.0)) / 120.0 * h)), time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * 5.0 / 12.0)), Scale(K3_iter, h * 5.0 / 12.0)), Scale(K4_iter, h / 12.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, (11.0 + sqrt(5.0)) / 120.0 * h)), Scale(K2_iter, (25.0 - sqrt(5.0)) / 120.0 * h)), Scale(K3_iter, (25.0 - 13.0 * sqrt(5.0)) / 120.0 * h)), Scale(K4_iter, (-1.0 + sqrt(5.0)) / 120.0 * h)), time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, (11.0 - sqrt(5.0)) / 120.0 * h)), Scale(K2_iter, (25.0 + 13.0 * sqrt(5.0)) / 120.0 * h)), Scale(K3_iter, (25.0 + sqrt(5.0)) / 120.0 * h)), Scale(K4_iter, (-1.0 - sqrt(5.0)) / 120.0 * h)), time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * 5.0 / 12.0)), Scale(K3_iter, h * 5.0 / 12.0)), Scale(K4_iter, h / 12.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Scale(K1_iter_n, 1.0 / 12.0), Scale(K2_iter_n, 5.0 / 12.0)), Scale(K3_iter_n, 5.0 / 12.0)), Scale(K4_iter_n, 1.0 / 12.0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 5) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
			vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (119.0 + 3 * sqrt(21.0)) / 1960.0)), Scale(K2_iter, h * (343.0 - 9.0 * sqrt(21.0)) / 2520.0)), Scale(K3_iter, h * (392.0 - 96.0 * sqrt(21.0)) / 2205.0)), Scale(K4_iter, h * (343.0 - 69.0 * sqrt(21.0)) / 2520.0)), Scale(K5_iter, h * (-21.0 + 3 * sqrt(21.0)) / 1960.0)), time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 13.0 / 320.0)), Scale(K2_iter, h * (392.0 + 105.0 * sqrt(21.0)) / 2880.0)), Scale(K3_iter, h * 8.0 / 45.0)), Scale(K4_iter, h * (392.0 - 105.0 * sqrt(21.0)) / 2880.0)), Scale(K5_iter, h * 3.0 / 320.0)), time[i - 1] + 0.5 * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (119.0 - 3 * sqrt(21.0)) / 1960.0)), Scale(K2_iter, h * (343.0 + 69.0 * sqrt(21.0)) / 2520.0)), Scale(K3_iter, h * (392.0 + 96.0 * sqrt(21.0)) / 2205.0)), Scale(K4_iter, h * (343.0 + 9.0 * sqrt(21.0)) / 2520.0)), Scale(K5_iter, h * (-21.0 - 3 * sqrt(21.0)) / 1960.0)), time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
			vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * 49.0 / 180.0)), Scale(K3_iter, h * 16.0 / 45.0)), Scale(K4_iter, h * 49.0 / 180.0)), Scale(K5_iter, h /20.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10 || Tolerance(K5_iter, K5_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K5_iter = K5_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (119.0 + 3 * sqrt(21.0)) / 1960.0)), Scale(K2_iter, h * (343.0 - 9.0 * sqrt(21.0)) / 2520.0)), Scale(K3_iter, h * (392.0 - 96.0 * sqrt(21.0)) / 2205.0)), Scale(K4_iter, h * (343.0 - 69.0 * sqrt(21.0)) / 2520.0)), Scale(K5_iter, h * (-21.0 + 3 * sqrt(21.0)) / 1960.0)), time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 13.0 / 320.0)), Scale(K2_iter, h * (392.0 + 105.0 * sqrt(21.0)) / 2880.0)), Scale(K3_iter, h * 8.0 / 45.0)), Scale(K4_iter, h * (392.0 - 105.0 * sqrt(21.0)) / 2880.0)), Scale(K5_iter, h * 3.0 / 320.0)), time[i - 1] + 0.5 * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (119.0 - 3 * sqrt(21.0)) / 1960.0)), Scale(K2_iter, h * (343.0 + 69.0 * sqrt(21.0)) / 2520.0)), Scale(K3_iter, h * (392.0 + 96.0 * sqrt(21.0)) / 2205.0)), Scale(K4_iter, h * (343.0 + 9.0 * sqrt(21.0)) / 2520.0)), Scale(K5_iter, h * (-21.0 - 3 * sqrt(21.0)) / 1960.0)), time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
				K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * 49.0 / 180.0)), Scale(K3_iter, h * 16.0 / 45.0)), Scale(K4_iter, h * 49.0 / 180.0)), Scale(K5_iter, h / 20.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Add(Scale(K1_iter_n, 1.0 / 20.0), Scale(K2_iter_n, 49.0 / 180.0)), Scale(K3_iter_n, 16.0 / 45.0)), Scale(K4_iter_n, 49.0 / 180.0)), Scale(K5_iter_n, 1.0 / 20.0)), h), x_solved[i - 1]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Lobatto IIIB  V
vector<vector <double>> ODE_Solve_RK_LIIIB(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	if (s == 2) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.5 * h)), Scale(K2_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 2.0)), Scale(K2_iter, 0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.5 * h)), Scale(K2_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 2.0)), Scale(K2_iter, 0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Add(Scale(K1_iter_n, h / 2.0), Scale(K2_iter_n, h / 2.0)), x_solved[i - 1]);
		}
	}
	else if (s == 3) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, - h / 6.0)), Scale(K3_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, h / 3.0)), Scale(K3_iter, 0)), time[i - 1] + 0.5 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 5.0 * h / 6.0)), Scale(K3_iter, 0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, -h / 6.0)), Scale(K3_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, h / 3.0)), Scale(K3_iter, 0)), time[i - 1] + 0.5 * h, K);
				K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 5.0 * h / 6.0)), Scale(K3_iter, 0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Scale(K1_iter_n, 1.0 / 6.0), Scale(K2_iter_n, 2.0 / 3.0)), Scale(K3_iter_n, 1.0 / 6.0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 4) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * (-1.0 - sqrt(5.0)) / 24.0)), Scale(K3_iter, h * (-1.0 + sqrt(5.0)) / 24.0)), Scale(K4_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, (25.0 + sqrt(5.0)) / 120.0 * h)), Scale(K3_iter, (25.0 - 13.0 * sqrt(5.0)) / 120.0 * h)), Scale(K4_iter, 0)), time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, (25.0 + 13.0 * sqrt(5.0)) / 120.0 * h)), Scale(K3_iter, (25.0 - sqrt(5.0)) / 120.0 * h)), Scale(K4_iter, 0)), time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * (11.0 - sqrt(5.0)) / 24.0)), Scale(K3_iter, h * (11.0 + sqrt(5.0)) / 24.0)), Scale(K4_iter, 0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * (-1.0 - sqrt(5.0)) / 24.0)), Scale(K3_iter, h * (-1.0 + sqrt(5.0)) / 24.0)), Scale(K4_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, (25.0 + sqrt(5.0)) / 120.0 * h)), Scale(K3_iter, (25.0 - 13.0 * sqrt(5.0)) / 120.0 * h)), Scale(K4_iter, 0)), time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, (25.0 + 13.0 * sqrt(5.0)) / 120.0 * h)), Scale(K3_iter, (25.0 - sqrt(5.0)) / 120.0 * h)), Scale(K4_iter, 0)), time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * (11.0 - sqrt(5.0)) / 24.0)), Scale(K3_iter, h * (11.0 + sqrt(5.0)) / 24.0)), Scale(K4_iter, 0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Scale(K1_iter_n, 1.0 / 12.0), Scale(K2_iter_n, 5.0 / 12.0)), Scale(K3_iter_n, 5.0 / 12.0)), Scale(K4_iter_n, 1.0 / 12.0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 5) { //doubt
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
			vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h* (-7.0 - sqrt(21.0)) / 120.0)), Scale(K3_iter, h / 15.0)), Scale(K4_iter, (-7.0 + sqrt(21.0)) / 120.0)), Scale(K5_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (343.0 + 9.0 * sqrt(21.0)) / 2520.0)), Scale(K3_iter, h * (56.0 - 15.0 * sqrt(21.0)) / 315.0)), Scale(K4_iter, h * (343.0 - 69.0 * sqrt(21.0)) / 2520.0)), Scale(K5_iter, 0)), time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (49.0 + 12.0 * sqrt(12.0)) / 360.0)), Scale(K3_iter, h * 8.0 / 45.0)), Scale(K4_iter, h * (49.0 - 12.0 * sqrt(12.0)) / 360.0)), Scale(K5_iter, 0)), time[i - 1] + 0.5 * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (343.0 + 69.0 * sqrt(21.0)) / 2520.0)), Scale(K3_iter, h * (56.0 + 15.0 * sqrt(21.0)) / 315.0)), Scale(K4_iter, h * (343.0 - 9.0 * sqrt(21.0)) / 2520.0)), Scale(K5_iter, 0)), time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
			vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (119.0 - 3.0 * sqrt(21.0)) / 360.0)), Scale(K3_iter, h * 13.0 / 45.0)), Scale(K4_iter, h * (119.0 + 3.0 * sqrt(21.0)) / 360.0)), Scale(K5_iter, 0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10 || Tolerance(K5_iter, K5_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K5_iter = K5_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (-7.0 - sqrt(21.0)) / 120.0)), Scale(K3_iter, h / 15.0)), Scale(K4_iter, (-7.0 + sqrt(21.0)) / 120.0)), Scale(K5_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (343.0 + 9.0 * sqrt(21.0)) / 2520.0)), Scale(K3_iter, h * (56.0 - 15.0 * sqrt(21.0)) / 315.0)), Scale(K4_iter, h * (343.0 - 69.0 * sqrt(21.0)) / 2520.0)), Scale(K5_iter, 0)), time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (49.0 + 12.0 * sqrt(12.0)) / 360.0)), Scale(K3_iter, h * 8.0 / 45.0)), Scale(K4_iter, h * (49.0 - 12.0 * sqrt(12.0)) / 360.0)), Scale(K5_iter, 0)), time[i - 1] + 0.5 * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (343.0 + 69.0 * sqrt(21.0)) / 2520.0)), Scale(K3_iter, h * (56.0 + 15.0 * sqrt(21.0)) / 315.0)), Scale(K4_iter, h * (343.0 - 9.0 * sqrt(21.0)) / 2520.0)), Scale(K5_iter, 0)), time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
				K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (119.0 - 3.0 * sqrt(21.0)) / 360.0)), Scale(K3_iter, h * 13.0 / 45.0)), Scale(K4_iter, h * (119.0 + 3.0 * sqrt(21.0)) / 360.0)), Scale(K5_iter, 0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Add(Scale(K1_iter_n, 1.0 / 20.0), Scale(K2_iter_n, 49.0 / 180.0)), Scale(K3_iter_n, 16.0 / 45.0)), Scale(K4_iter_n, 49.0 / 180.0)), Scale(K5_iter_n, 1.0 / 20.0)), h), x_solved[i - 1]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Lobatto IIIC  V
vector<vector <double>> ODE_Solve_RK_LIIIC(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	if (s == 2) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.5 * h)), Scale(K2_iter, -0.5 * h)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 2.0)), Scale(K2_iter, 0.5 * h)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.5 * h)), Scale(K2_iter, -0.5 * h)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 2.0)), Scale(K2_iter, 0.5 * h)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Add(Scale(K1_iter_n, h / 2.0), Scale(K2_iter_n, h / 2.0)), x_solved[i - 1]);
		}
	}
	else if (s == 3) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, -h / 3.0)), Scale(K3_iter, h / 6.0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, h * 5.0 / 12.0)), Scale(K3_iter, - h /12.0)), time[i - 1] + 0.5 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 2.0 * h / 3.0)), Scale(K3_iter, h / 6.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, -h / 3.0)), Scale(K3_iter, h / 6.0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, h * 5.0 / 12.0)), Scale(K3_iter, -h / 12.0)), time[i - 1] + 0.5 * h, K);
				K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 2.0 * h / 3.0)), Scale(K3_iter, h / 6.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Scale(K1_iter_n, 1.0 / 6.0), Scale(K2_iter_n, 2.0 / 3.0)), Scale(K3_iter_n, 1.0 / 6.0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 4) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * -sqrt(5.0) / 12.0)), Scale(K3_iter, h * sqrt(5.0) / 12.0)), Scale(K4_iter, -h / 12.0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h / 4.0)), Scale(K3_iter, (10.0 - 7.0 * sqrt(5.0)) / 60.0 * h)), Scale(K4_iter, h * sqrt(5.0) / 60.0)), time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, (10.0 + 7.0 * sqrt(5.0)) / 60.0 * h)), Scale(K3_iter, h / 4.0)), Scale(K4_iter, h * -sqrt(5.0) / 60.0)), time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * 5.0 / 12.0)), Scale(K3_iter, h * 5.0 / 12.0)), Scale(K4_iter, h /12.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * -sqrt(5.0) / 12.0)), Scale(K3_iter, h * sqrt(5.0) / 12.0)), Scale(K4_iter, -h / 12.0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h / 4.0)), Scale(K3_iter, (10.0 - 7.0 * sqrt(5.0)) / 60.0 * h)), Scale(K4_iter, h * sqrt(5.0) / 60.0)), time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, (10.0 + 7.0 * sqrt(5.0)) / 60.0 * h)), Scale(K3_iter, h / 4.0)), Scale(K4_iter, h * -sqrt(5.0) / 60.0)), time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 12.0)), Scale(K2_iter, h * 5.0 / 12.0)), Scale(K3_iter, h * 5.0 / 12.0)), Scale(K4_iter, h / 12.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Scale(K1_iter_n, 1.0 / 12.0), Scale(K2_iter_n, 5.0 / 12.0)), Scale(K3_iter_n, 5.0 / 12.0)), Scale(K4_iter_n, 1.0 / 12.0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 5) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
			vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * - 7.0 / 60.0)), Scale(K3_iter, h * 2.0 / 15.0)), Scale(K4_iter, h * -7.0  / 60.0)), Scale(K5_iter, h / 20.0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * 29.0 / 180.0)), Scale(K3_iter, h * (47.0 - 15.0 * sqrt(21.0)) / 315.0)), Scale(K4_iter, h * (203.0 - 30.0 * sqrt(21.0)) / 1260.0)), Scale(K5_iter, - 3.0 / 140.0 * h)), time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (329.0 + 105.0 * sqrt(21.0)) / 2880.0)), Scale(K3_iter, h * 73.0 / 360.0)), Scale(K4_iter, h * (329.0 - 105.0 * sqrt(21.0)) / 2880.0)), Scale(K5_iter, 3.0 / 160.0 * h)), time[i - 1] + 0.5 * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (203.0 + 30.0 * sqrt(21.0)) / 1260.0)), Scale(K3_iter, h * (47.0 + 15.0 * sqrt(21.0)) / 315.0)), Scale(K4_iter, h * 29.0 / 180.0)), Scale(K5_iter, h * -3.0/140.0)), time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
			vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * 49.0 / 180.0)), Scale(K3_iter, h * 16.0 / 45.0)), Scale(K4_iter, h * 49.0 / 180.0)), Scale(K5_iter, h / 20.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10 || Tolerance(K5_iter, K5_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K5_iter = K5_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * -7.0 / 60.0)), Scale(K3_iter, h * 2.0 / 15.0)), Scale(K4_iter, h * -7.0 / 60.0)), Scale(K5_iter, h / 20.0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * 29.0 / 180.0)), Scale(K3_iter, h * (47.0 - 15.0 * sqrt(21.0)) / 315.0)), Scale(K4_iter, h * (203.0 - 30.0 * sqrt(21.0)) / 1260.0)), Scale(K5_iter, -3.0 / 140.0 * h)), time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (329.0 + 105.0 * sqrt(21.0)) / 2880.0)), Scale(K3_iter, h * 73.0 / 360.0)), Scale(K4_iter, h * (329.0 - 105.0 * sqrt(21.0)) / 2880.0)), Scale(K5_iter, 3.0 / 160.0 * h)), time[i - 1] + 0.5 * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * (203.0 + 30.0 * sqrt(21.0)) / 1260.0)), Scale(K3_iter, h * (47.0 + 15.0 * sqrt(21.0)) / 315.0)), Scale(K4_iter, h * 29.0 / 180.0)), Scale(K5_iter, h * -3.0 / 140.0)), time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
				K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 20.0)), Scale(K2_iter, h * 49.0 / 180.0)), Scale(K3_iter, h * 16.0 / 45.0)), Scale(K4_iter, h * 49.0 / 180.0)), Scale(K5_iter, h / 20.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Add(Scale(K1_iter_n, 1.0 / 20.0), Scale(K2_iter_n, 49.0 / 180.0)), Scale(K3_iter_n, 16.0 / 45.0)), Scale(K4_iter_n, 49.0 / 180.0)), Scale(K5_iter_n, 1.0 / 20.0)), h), x_solved[i - 1]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Lobatto IIIC*  V
vector<vector <double>> ODE_Solve_RK_LIIICn(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	if (s == 2) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h)), Scale(K2_iter, 0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h)), Scale(K2_iter, 0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Add(Scale(K1_iter_n, h / 2.0), Scale(K2_iter_n, h / 2.0)), x_solved[i - 1]);
		}
	}
	else if (s == 3) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 4.0)), Scale(K2_iter, h / 4.0)), Scale(K3_iter, 0)), time[i - 1] + 0.5 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, h)), Scale(K3_iter, 0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 4.0)), Scale(K2_iter, h / 4.0)), Scale(K3_iter, 0)), time[i - 1] + 0.5 * h, K);
				K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, h)), Scale(K3_iter, 0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Scale(K1_iter_n, 1.0 / 6.0), Scale(K2_iter_n, 2.0 / 3.0)), Scale(K3_iter_n, 1.0 / 6.0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 4) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (5.0 + sqrt(5.0)) / 60.0)), Scale(K2_iter, h / 6.0)), Scale(K3_iter, (15.0 - 7.0 * sqrt(5.0)) / 60.0 * h)), Scale(K4_iter, 0)), time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (5.0 - sqrt(5.0)) / 60.0)), Scale(K2_iter, (15.0 + 7.0 * sqrt(5.0)) / 60.0 * h)), Scale(K3_iter, h / 6.0)), Scale(K4_iter, 0)), time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, h * (5.0-sqrt(5.0)) / 12.0)), Scale(K3_iter, h * (5.0 + sqrt(5.0)) / 12.0)), Scale(K4_iter, 0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (5.0 + sqrt(5.0)) / 60.0)), Scale(K2_iter, h / 6.0)), Scale(K3_iter, (15.0 - 7.0 * sqrt(5.0)) / 60.0 * h)), Scale(K4_iter, 0)), time[i - 1] + (0.5 - sqrt(5.0) / 10.0) * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (5.0 - sqrt(5.0)) / 60.0)), Scale(K2_iter, (15.0 + 7.0 * sqrt(5.0)) / 60.0 * h)), Scale(K3_iter, h / 6.0)), Scale(K4_iter, 0)), time[i - 1] + (0.5 + sqrt(5.0) / 10.0) * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, h * (5.0 - sqrt(5.0)) / 12.0)), Scale(K3_iter, h * (5.0 + sqrt(5.0)) / 12.0)), Scale(K4_iter, 0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Scale(K1_iter_n, 1.0 / 12.0), Scale(K2_iter_n, 5.0 / 12.0)), Scale(K3_iter_n, 5.0 / 12.0)), Scale(K4_iter_n, 1.0 / 12.0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 5) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
			vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 14.0)), Scale(K2_iter, h / 9.0)), Scale(K3_iter, h * (13.0 - 3.0 * sqrt(21.0)) / 63.0)), Scale(K4_iter, h * (14.0 - 3.0 * sqrt(21.0)) / 126.0)), Scale(K5_iter, 0)), time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 32.0)), Scale(K2_iter, h * (91.0 + 21.0 * sqrt(21.0)) / 576.0)), Scale(K3_iter, h * 11.0 / 72.0)), Scale(K4_iter, h * (91.0 - 21.0 * sqrt(21.0)) / 576.0)), Scale(K5_iter, 0)), time[i - 1] + 0.5 * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 14.0)), Scale(K2_iter, h * (14.0 + 3.0 * sqrt(21.0)) / 126.0)), Scale(K3_iter, h * (13.0 + 3.0 * sqrt(21.0)) / 63.0)), Scale(K4_iter, h / 9.0)), Scale(K5_iter, 0)), time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
			vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, h * 7.0 / 18.0)), Scale(K3_iter, h * 2.0 / 9.0)), Scale(K4_iter, h * 7.0 / 18.0)), Scale(K5_iter, 0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10 || Tolerance(K5_iter, K5_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K5_iter = K5_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 14.0)), Scale(K2_iter, h / 9.0)), Scale(K3_iter, h * (13.0 - 3.0 * sqrt(21.0)) / 63.0)), Scale(K4_iter, h * (14.0 - 3.0 * sqrt(21.0)) / 126.0)), Scale(K5_iter, 0)), time[i - 1] + (0.5 - sqrt(21.0) / 14.0) * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 32.0)), Scale(K2_iter, h * (91.0 + 21.0 * sqrt(21.0)) / 576.0)), Scale(K3_iter, h * 11.0 / 72.0)), Scale(K4_iter, h * (91.0 - 21.0 * sqrt(21.0)) / 576.0)), Scale(K5_iter, 0)), time[i - 1] + 0.5 * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 14.0)), Scale(K2_iter, h * (14.0 + 3.0 * sqrt(21.0)) / 126.0)), Scale(K3_iter, h * (13.0 + 3.0 * sqrt(21.0)) / 63.0)), Scale(K4_iter, h / 9.0)), Scale(K5_iter, 0)), time[i - 1] + (0.5 + sqrt(21.0) / 14.0) * h, K);
				K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, h * 7.0 / 18.0)), Scale(K3_iter, h * 2.0 / 9.0)), Scale(K4_iter, h * 7.0 / 18.0)), Scale(K5_iter, 0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Add(Scale(K1_iter_n, 1.0 / 20.0), Scale(K2_iter_n, 49.0 / 180.0)), Scale(K3_iter_n, 16.0 / 45.0)), Scale(K4_iter_n, 49.0 / 180.0)), Scale(K5_iter_n, 1.0 / 20.0)), h), x_solved[i - 1]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K SDIRK  V
vector<vector <double>> ODE_Solve_RK_SDIRK(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, int s, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	if (s == 2) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1] + h * (1.0 - sqrt(2.0) / 2.0), K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (1.0 - sqrt(2.0) / 2.0))), Scale(K2_iter, 0)), time[i - 1] + h * (1.0 - sqrt(2.0) / 2.0), K);
			vector<double> K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h * sqrt(2.0) / 2.0)), Scale(K2_iter, h * (1.0 - sqrt(2.0) / 2.0))), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K1_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (1.0 - sqrt(2.0) / 2.0))), Scale(K2_iter, 0)), time[i - 1] + h * (1.0 - sqrt(2.0) / 2.0), K);
				K2_iter_n = Fun(Add(Add(x_solved[i - 1], Scale(K1_iter, h * sqrt(2.0) / 2.0)), Scale(K2_iter, h * (1.0 - sqrt(2.0) / 2.0))), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Add(Scale(K1_iter_n, h * sqrt(2.0) / 2.0), Scale(K2_iter_n, h * (1.0 - sqrt(2.0) / 2.0))), x_solved[i - 1]);
		}
	}
	else if (s == 3) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1] + h * (3.0 - sqrt(3.0)) / 6.0, K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + h * (1.0 - sqrt(3.0) / 3.0), K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (3.0 - sqrt(3.0)) / 6.0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (3.0 - sqrt(3.0)) / 6.0)), Scale(K2_iter, h * (3.0 - sqrt(3.0)) / 6.0)), Scale(K3_iter, 0)), time[i - 1] + 0.5 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (3.0 - sqrt(3.0)) / 6.0)), Scale(K2_iter, h * sqrt(3.0) / 3.0)), Scale(K3_iter, h * (3.0 - sqrt(3.0)) / 6.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K1_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (3.0 - sqrt(3.0)) / 6.0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (3.0 - sqrt(3.0)) / 6.0)), Scale(K2_iter, h * (3.0 - sqrt(3.0)) / 6.0)), Scale(K3_iter, 0)), time[i - 1] + 0.5 * h, K);
				K3_iter_n = Fun(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (3.0 - sqrt(3.0)) / 6.0)), Scale(K2_iter, h * sqrt(3.0) / 3.0)), Scale(K3_iter, h * (3.0 - sqrt(3.0)) / 6.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Scale(K1_iter_n, (3.0 - sqrt(3.0)) / 6.0), Scale(K2_iter_n, sqrt(3.0) / 3.0)), Scale(K3_iter_n, (3.0 - sqrt(3.0)) / 6.0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 4) {// doubt
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			
			double c1 = (3.0 - sqrt(3.0)) / 6.0;
			double thi = 0.5 + 2.0 * sqrt(3.0) / 9.0;

			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + c1 * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + (1.0 - c1) * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);
			
			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), time[i - 1], K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (6.0 * (c1 + thi) - 5.0) / 12.0)), Scale(K2_iter, h * (1.0 - thi) / 2.0)), Scale(K3_iter, h * (1.0 - thi) / 2.0)), Scale(K4_iter, (6.0 * (c1 + thi) - 7.0) / 12.0)), time[i - 1] + c1 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (7.0 - 6.0*(c1 + thi)) / 12.0)), Scale(K2_iter, h * thi / 2.0)), Scale(K3_iter, h * thi / 2.0)), Scale(K4_iter, (5.0 - 6.0 * (c1+thi)) / 12.0)), time[i - 1] + (1.0 - c1) * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, h / 2.0)), Scale(K3_iter, h / 2.0)), Scale(K4_iter, 0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), time[i - 1], K);
				K2_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (6.0 * (c1 + thi) - 5.0) / 12.0)), Scale(K2_iter, h * (1.0 - thi) / 2.0)), Scale(K3_iter, h * (1.0 - thi) / 2.0)), Scale(K4_iter, (6.0 * (c1 + thi) - 7.0) / 12.0)), time[i - 1] + c1 * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * (7.0 - 6.0 * (c1 + thi)) / 12.0)), Scale(K2_iter, h * thi / 2.0)), Scale(K3_iter, h * thi / 2.0)), Scale(K4_iter, (5.0 - 6.0 * (c1 + thi)) / 12.0)), time[i - 1] + (1.0 - c1) * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, h / 2.0)), Scale(K3_iter, h / 2.0)), Scale(K4_iter, 0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Scale(K1_iter_n, 0), Scale(K2_iter_n, 0.5)), Scale(K3_iter_n, 0.5)), Scale(K4_iter_n, 0)), h), x_solved[i - 1]);
		}
	}
	else if (s == 5) {
		for (int i = 1; i < len; ++i) {
			time[i] = time[i - 1] + h;
			vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1] + 0.25 * h, K);
			vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 0.75 * h, K);
			vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + 11.0 / 20.0 * h, K);
			vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
			vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

			vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.25 * h)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1] + 0.25 * h, K);
			vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 2.0)), Scale(K2_iter, h / 4.0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1] + 0.75 * h, K);
			vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 17.0 / 50.0)), Scale(K2_iter, - h / 25.0)), Scale(K3_iter, h  / 4.0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1] + 11.0 / 20.0 * h, K);
			vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 371.0 / 1360.0)), Scale(K2_iter, h * - 137.0 / 2720.0)), Scale(K3_iter, h * 15.0 / 544.0)), Scale(K4_iter, h / 4.0)), Scale(K5_iter, 0)), time[i - 1] + 0.5 * h, K);
			vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 25.0 / 24.0)), Scale(K2_iter, h * - 49.0 / 48.0)), Scale(K3_iter, h * 125.0 / 16.0)), Scale(K4_iter, h * - 85.0 / 12.0)), Scale(K5_iter, h / 4.0)), time[i - 1] + h, K);

			while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10 || Tolerance(K5_iter, K5_iter_n) > 1e-10) {
				K1_iter = K1_iter_n;
				K2_iter = K2_iter_n;
				K3_iter = K3_iter_n;
				K4_iter = K4_iter_n;
				K5_iter = K5_iter_n;
				K1_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.25 * h)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1] + 0.25 * h, K);
				K2_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 2.0)), Scale(K2_iter, h / 4.0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1] + 0.75 * h, K);
				K3_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 17.0 / 50.0)), Scale(K2_iter, -h / 25.0)), Scale(K3_iter, h / 4.0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), time[i - 1] + 11.0 / 20.0 * h, K);
				K4_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 371.0 / 1360.0)), Scale(K2_iter, h * -137.0 / 2720.0)), Scale(K3_iter, h * 15.0 / 544.0)), Scale(K4_iter, h / 4.0)), Scale(K5_iter, 0)), time[i - 1] + 0.5 * h, K);
				K5_iter_n = Fun(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 25.0 / 24.0)), Scale(K2_iter, h * -49.0 / 48.0)), Scale(K3_iter, h * 125.0 / 16.0)), Scale(K4_iter, h * -85.0 / 12.0)), Scale(K5_iter, h / 4.0)), time[i - 1] + h, K);
			}

			x_solved[i] = Add(Scale(Add(Add(Add(Add(Scale(K1_iter_n, 25.0 / 24.0), Scale(K2_iter_n, - 49.0 / 48.0)), Scale(K3_iter_n, 125.0 / 16.0)), Scale(K4_iter_n, - 85.0 / 12.0)), Scale(K5_iter_n, 1.0 / 4.0)), h), x_solved[i - 1]);
		}
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K ESDIRK  Exp.......
vector<vector <double>> ODE_Solve_RK_ESDIRK(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;
		
	for (int i = 1; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 0.5 * h, K);
		vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + (2.0 + sqrt(2.0)) / 4.0 * h, K);
		vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + 0.53 * h, K);
		vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + 0.8 * h, K);
		vector<double> K6_iter = Fun(x_solved[i - 1], time[i - 1] + 17.0 / 25.0 * h, K);
		vector<double> K7_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);
		vector<double> K8_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

		vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1], K);
		vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.25 * h)), Scale(K2_iter, 0.25 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.5 * h, K);
		vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1748874742213.0 / 5795261096931.0 * h)), Scale(K2_iter, 1748874742213.0 / 5795261096931.0 * h)), Scale(K3_iter, 0.25 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + (2.0 + sqrt(2.0)) / 4.0 * h, K);
		vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 2426486750897.0 / 12677310711630.0 * h)), Scale(K2_iter, 2426486750897.0 / 12677310711630.0 * h)), Scale(K3_iter, -783385356511.0 / 7619901499812.0 * h)), Scale(K4_iter, 0.25 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.53 * h, K);
		vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1616209367427.0 / 5722977998639.0 * h)), Scale(K2_iter, 1616209367427.0 / 5722977998639.0 * h)), Scale(K3_iter, -211896077633.0 / 5134769641545.0 * h)), Scale(K4_iter, 464248917192.0 / 17550087120101.0 * h)), Scale(K5_iter, 0.25 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.8 * h, K);
		vector<double> K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1860464898611.0 / 7805430689312.0 * h)), Scale(K2_iter, 1825204367749.0 / 7149715425471.0 * h)), Scale(K3_iter, -1289376786583.0 / 6598860380111.0 * h)), Scale(K4_iter, 55566826943.0 / 2961051076052.0 * h)), Scale(K5_iter, 1548994872005.0 / 13709222415197.0 * h)), Scale(K6_iter, 0.25 * h)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 17.0 / 25.0 * h, K);
		vector<double> K7_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1783640092711.0 / 14417713428467.0 * h)), Scale(K2_iter, -5781183663275.0 / 18946039887294.0 * h)), Scale(K3_iter, 57847255876685.0 / 10564937217081.0 * h)), Scale(K4_iter, 29339178902168.0 / 9787613280015.0 * h)), Scale(K5_iter, 122011506936853.0 / 12523522131766.0 * h)), Scale(K6_iter, -60418758964762.0 / 9539790648093.0 * h)), Scale(K7_iter, 0.25 * h)), Scale(K8_iter, 0)), time[i - 1] + h, K);
		vector<double> K8_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 3148564786223.0 / 23549948766475.0 * h)), Scale(K2_iter, -4152366519273.0 / 20368318839251.0 * h)), Scale(K3_iter, -143958253112335.0 / 33767350176582.0 * h)), Scale(K4_iter, 16929685656751.0 / 6821330976083.0 * h)), Scale(K5_iter, 37330861322165.0 / 4907624269821.0 * h)), Scale(K6_iter, -103974720808012.0 / 20856851060343.0 * h)), Scale(K7_iter, -93596557767.0 / 4675692258479.0 * h)), Scale(K8_iter, 0.25 * h)), time[i - 1] + h, K);

		while (Tolerance(K1_iter, K1_iter_n) > 1e-20 || Tolerance(K2_iter, K2_iter_n) > 1e-20 || Tolerance(K3_iter, K3_iter_n) > 1e-20 || Tolerance(K4_iter, K4_iter_n) > 1e-20 || Tolerance(K5_iter, K5_iter_n) > 1e-20 || Tolerance(K6_iter, K6_iter_n) > 1e-20 || Tolerance(K7_iter, K7_iter_n) > 1e-20 || Tolerance(K8_iter, K8_iter_n) > 1e-20) {
			K1_iter = K1_iter_n;
			K2_iter = K2_iter_n;
			K3_iter = K3_iter_n;
			K4_iter = K4_iter_n;
			K5_iter = K5_iter_n;
			K6_iter = K6_iter_n;
			K7_iter = K7_iter_n;
			K8_iter = K8_iter_n;

			K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1], K);
			K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.25 * h)), Scale(K2_iter, 0.25 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.5 * h, K);
			K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1748874742213.0 / 5795261096931.0 * h)), Scale(K2_iter, 1748874742213.0 / 5795261096931.0 * h)), Scale(K3_iter, 0.25 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + (2.0 + sqrt(2.0)) / 4.0 * h, K);
			K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 2426486750897.0 / 12677310711630.0 * h)), Scale(K2_iter, 2426486750897.0 / 12677310711630.0 * h)), Scale(K3_iter, -783385356511.0 / 7619901499812.0 * h)), Scale(K4_iter, 0.25 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.53 * h, K);
			K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1616209367427.0 / 5722977998639.0 * h)), Scale(K2_iter, 1616209367427.0 / 5722977998639.0 * h)), Scale(K3_iter, -211896077633.0 / 5134769641545.0 * h)), Scale(K4_iter, 464248917192.0 / 17550087120101.0 * h)), Scale(K5_iter, 0.25 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.8 * h, K);
			K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1860464898611.0 / 7805430689312.0 * h)), Scale(K2_iter, 1825204367749.0 / 7149715425471.0 * h)), Scale(K3_iter, -1289376786583.0 / 6598860380111.0 * h)), Scale(K4_iter, 55566826943.0 / 2961051076052.0 * h)), Scale(K5_iter, 1548994872005.0 / 13709222415197.0 * h)), Scale(K6_iter, 0.25 * h)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 17.0 / 25.0 * h, K);
			K7_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1783640092711.0 / 14417713428467.0 * h)), Scale(K2_iter, -5781183663275.0 / 18946039887294.0 * h)), Scale(K3_iter, 57847255876685.0 / 10564937217081.0 * h)), Scale(K4_iter, 29339178902168.0 / 9787613280015.0 * h)), Scale(K5_iter, 122011506936853.0 / 12523522131766.0 * h)), Scale(K6_iter, -60418758964762.0 / 9539790648093.0 * h)), Scale(K7_iter, 0.25 * h)), Scale(K8_iter, 0)), time[i - 1] + h, K);
			K8_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 3148564786223.0 / 23549948766475.0 * h)), Scale(K2_iter, -4152366519273.0 / 20368318839251.0 * h)), Scale(K3_iter, -143958253112335.0 / 33767350176582.0 * h)), Scale(K4_iter, 16929685656751.0 / 6821330976083.0 * h)), Scale(K5_iter, 37330861322165.0 / 4907624269821.0 * h)), Scale(K6_iter, -103974720808012.0 / 20856851060343.0 * h)), Scale(K7_iter, -93596557767.0 / 4675692258479.0 * h)), Scale(K8_iter, 0.25 * h)), time[i - 1] + h, K);
		}

		x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Scale(K1_iter_n, 1707112744407.0 / 11125291145125.0), Scale(K2_iter_n, -34114578494.0 / 9511465441463.0)), Scale(K3_iter_n, -10730340352595.0 / 5750059075211.0)), Scale(K4_iter_n, 16308974155447.0 / 11154981868028.0)), Scale(K5_iter_n, 16015983083570.0 / 4734398780449.0)), Scale(K6_iter_n, -16745095336747.0 / 7220642436550.0)), Scale(K7_iter_n, -3941055932791.0 / 7113931146175.0)), Scale(K8_iter_n, 9504350599766.0 / 12768012822991.0)), h), x_solved[i - 1]);
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Khashin  Exp..........
vector<vector <double>> ODE_Solve_RKK(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 0.4 * h, K);
		vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + 0.6 * h, K);
		vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + h / 7.0, K);
		vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + h / 3.0, K);
		vector<double> K6_iter = Fun(x_solved[i - 1], time[i - 1] + 17.0 / 23.0 * h, K);
		vector<double> K7_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);
		vector<double> K8_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

		vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1], K);
		vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.4 * h)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.4 * h, K);
		vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 7.0 / 20.0 * h)), Scale(K2_iter, 1.0 / 4.0 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.6 * h, K);
		vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, - 659.0 / 1372.0 * h)), Scale(K2_iter, 1431.0 / 1372.0 * h)), Scale(K3_iter, -144.0 / 343.0 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h / 7.0, K);
		vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 91.0 / 324.0 * h)), Scale(K2_iter, 73.0 / 2916.0 * h)), Scale(K3_iter, -1.0 / 162.0 * h)), Scale(K4_iter, 49.0 / 1458.0 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h / 3.0, K);
		vector<double> K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -259.0 / 48668.0 * h)), Scale(K2_iter, -10653149.0 / 3358092.0 * h)), Scale(K3_iter, 217152.0 / 279841.0 * h)), Scale(K4_iter, -193648.0 / 839523.0 * h)), Scale(K5_iter, 943488.0 / 279841.0 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 17.0 / 23.0 * h, K);
		vector<double> K7_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 4885.0 / 702.0 * h)), Scale(K3_iter, -2855.0 / 1536.0 * h)), Scale(K4_iter, 22099.0 / 27648.0 * h)), Scale(K5_iter, -1359.0 / 224.0 * h)), Scale(K6_iter, 279841.0 / 279552.0 * h)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h, K);
		vector<double> K8_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -5.0 / 4.0 * h)), Scale(K2_iter, -1.0 / 3.0 * h)), Scale(K3_iter, - 25293.0 / 4096.0 * h)), Scale(K4_iter, -51401.0 / 122880.0 * h)), Scale(K5_iter, 43821.0 / 8960.0 * h)), Scale(K6_iter, -1228867.0 / 286720.0 * h)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h, K);

		while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10 || Tolerance(K5_iter, K5_iter_n) > 1e-10 || Tolerance(K6_iter, K6_iter_n) > 1e-10 || Tolerance(K7_iter, K7_iter_n) > 1e-10 || Tolerance(K8_iter, K8_iter_n) > 1e-10) {
			K1_iter = K1_iter_n;
			K2_iter = K2_iter_n;
			K3_iter = K3_iter_n;
			K4_iter = K4_iter_n;
			K5_iter = K5_iter_n;
			K6_iter = K6_iter_n;
			K7_iter = K7_iter_n;
			K8_iter = K8_iter_n;

			K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1], K);
			K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0.4 * h)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.4 * h, K);
			K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 7.0 / 20.0 * h)), Scale(K2_iter, 1.0 / 4.0 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 0.6 * h, K);
			K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -659.0 / 1372.0 * h)), Scale(K2_iter, 1431.0 / 1372.0 * h)), Scale(K3_iter, -144.0 / 343.0 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h / 7.0, K);
			K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 91.0 / 324.0 * h)), Scale(K2_iter, 73.0 / 2916.0 * h)), Scale(K3_iter, -1.0 / 162.0 * h)), Scale(K4_iter, 49.0 / 1458.0 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h / 3.0, K);
			K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -259.0 / 48668.0 * h)), Scale(K2_iter, -10653149.0 / 3358092.0 * h)), Scale(K3_iter, 217152.0 / 279841.0 * h)), Scale(K4_iter, -193648.0 / 839523.0 * h)), Scale(K5_iter, 943488.0 / 279841.0 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 17.0 / 23.0 * h, K);
			K7_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 4885.0 / 702.0 * h)), Scale(K3_iter, -2855.0 / 1536.0 * h)), Scale(K4_iter, 22099.0 / 27648.0 * h)), Scale(K5_iter, -1359.0 / 224.0 * h)), Scale(K6_iter, 279841.0 / 279552.0 * h)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h, K);
			K8_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -5.0 / 4.0 * h)), Scale(K2_iter, -1.0 / 3.0 * h)), Scale(K3_iter, -25293.0 / 4096.0 * h)), Scale(K4_iter, -51401.0 / 122880.0 * h)), Scale(K5_iter, 43821.0 / 8960.0 * h)), Scale(K6_iter, -1228867.0 / 286720.0 * h)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h, K);
		}

		x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Scale(K1_iter_n, 1.0 / 12.0), Scale(K2_iter_n, -2375.0 / 16848.0)), Scale(K3_iter_n, -2375.0 / 98304.0)), Scale(K4_iter_n, 271313.0 / 5308416.0)), Scale(K5_iter_n, 7695.0 / 14336.0)), Scale(K6_iter_n, 67441681.0 / 161021952.0)), Scale(K7_iter_n, 1.0 / 12.0)), Scale(K8_iter_n, -19.0 / 2304.0)), h), x_solved[i - 1]);
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Verner  V
vector<vector <double>> ODE_Solve_RKV(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + h / 6.0, K);
		vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + 4.0 / 15.0 * h, K);
		vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + 2.0 * h / 3.0, K);
		vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + h * 5.0 / 6.0, K);
		vector<double> K6_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);
		vector<double> K7_iter = Fun(x_solved[i - 1], time[i - 1] + h / 15.0, K);
		vector<double> K8_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

		vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1], K);
		vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h / 6.0, K);
		vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 4.0 / 75.0 * h)), Scale(K2_iter, 16.0 / 75.0 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 4.0 / 15.0 * h, K);
		vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 5.0 / 6.0 * h)), Scale(K2_iter, -8.0 / 3.0 * h)), Scale(K3_iter, 5.0 / 2.0 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 2.0 * h / 3.0, K);
		vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -165.0 / 64.0 * h)), Scale(K2_iter, 55.0 / 6.0 * h)), Scale(K3_iter, -425.0 / 64.0 * h)), Scale(K4_iter, 85.0 / 96.0 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h * 5.0 / 6.0, K);
		vector<double> K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 12.0 / 5.0 * h)), Scale(K2_iter, -8.0 * h)), Scale(K3_iter, 4015.0 / 612.0 * h)), Scale(K4_iter, -11.0 / 36.0 * h)), Scale(K5_iter, 86.0 / 225.0 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h, K);
		vector<double> K7_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 8263.0 / -15000.0)), Scale(K2_iter, 124.0 / 75.0 * h)), Scale(K3_iter, -643.0/680.0 * h)), Scale(K4_iter, -84.0 / 250.0 * h)), Scale(K5_iter, 2484.0 / 10625.0 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h / 15.0, K);
		vector<double> K8_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 3501.0 / 1720.0 * h)), Scale(K2_iter, -300.0 / 43.0 * h)), Scale(K3_iter, 297278.0 / 52632.0 * h)), Scale(K4_iter, -319.0 / 2322.0 * h)), Scale(K5_iter, 24068.0 / 84065.0 * h)), Scale(K6_iter, 3850.0 / 26703.0 * h)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h, K);

		while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10 || Tolerance(K5_iter, K5_iter_n) > 1e-10 || Tolerance(K6_iter, K6_iter_n) > 1e-10 || Tolerance(K7_iter, K7_iter_n) > 1e-10 || Tolerance(K8_iter, K8_iter_n) > 1e-10) {
			K1_iter = K1_iter_n;
			K2_iter = K2_iter_n;
			K3_iter = K3_iter_n;
			K4_iter = K4_iter_n;
			K5_iter = K5_iter_n;
			K6_iter = K6_iter_n;
			K7_iter = K7_iter_n;
			K8_iter = K8_iter_n;

			K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1], K);
			K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h / 6.0, K);
			K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 4.0 / 75.0 * h)), Scale(K2_iter, 16.0 / 75.0 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 4.0 / 15.0 * h, K);
			K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 5.0 / 6.0 * h)), Scale(K2_iter, -8.0 / 3.0 * h)), Scale(K3_iter, 5.0 / 2.0 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 2.0 * h / 3.0, K);
			K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -165.0 / 64.0 * h)), Scale(K2_iter, 55.0 / 6.0 * h)), Scale(K3_iter, -425.0 / 64.0 * h)), Scale(K4_iter, 85.0 / 96.0 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h * 5.0 / 6.0, K);
			K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 12.0 / 5.0 * h)), Scale(K2_iter, -8.0 * h)), Scale(K3_iter, 4015.0 / 612.0 * h)), Scale(K4_iter, -11.0 / 36.0 * h)), Scale(K5_iter, 86.0 / 225.0 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h, K);
			K7_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 8263.0 / -15000.0)), Scale(K2_iter, 124.0 / 75.0 * h)), Scale(K3_iter, -643.0 / 680.0 * h)), Scale(K4_iter, -84.0 / 250.0 * h)), Scale(K5_iter, 2484.0 / 10625.0 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h / 15.0, K);
			K8_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 3501.0 / 1720.0 * h)), Scale(K2_iter, -300.0 / 43.0 * h)), Scale(K3_iter, 297278.0 / 52632.0 * h)), Scale(K4_iter, -319.0 / 2322.0 * h)), Scale(K5_iter, 24068.0 / 84065.0 * h)), Scale(K6_iter, 3850.0 / 26703.0 * h)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h, K);
		}

		x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Scale(K1_iter_n, 3.0 / 40.0), Scale(K2_iter_n, 0)), Scale(K3_iter_n, 875.0 / 2244.0)), Scale(K4_iter_n, 23.0 / 72.0)), Scale(K5_iter_n, 264.0 / 1955.0)), Scale(K6_iter_n, 0)), Scale(K7_iter_n, 125.0 / 11592.0)), Scale(K8_iter_n, 43.0 / 616.0)), h), x_solved[i - 1]);
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Tsitouras 5  V
vector<vector <double>> ODE_Solve_RKT5(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1] + 2.0 / 3.0 * h, K);
		vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 1.0 / 6.0 * h, K);
		vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + 0.75 * h, K);
		vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);
		vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + h * 0.8, K);
		vector<double> K6_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);
	
		vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + 2.0 / 3.0 * h, K);
		vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + 1.0 / 6.0 * h, K);
		vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, - 15.0 / 8.0 * h)), Scale(K2_iter, 21.0 / 8.0 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + 0.75 * h, K);
		vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -9.0 * h)), Scale(K2_iter, 75.0 / 7.0 * h)), Scale(K3_iter, -5.0 / 7.0 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + h, K);
		vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -3.0 * h)), Scale(K2_iter, 34257.0 / 8750.0 * h)), Scale(K3_iter, -114.0 / 875.0 * h)), Scale(K4_iter, 19.0 / 1250.0 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + h * 0.8, K);
		vector<double> K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 123.0 / 380.0 * h)), Scale(K3_iter, 2.5 * h)), Scale(K4_iter, 3.0 / 20.0 * h)), Scale(K5_iter, - 75.0 / 38.0 * h)), Scale(K6_iter, 0)), time[i - 1] + h, K);
		
		while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10 || Tolerance(K5_iter, K5_iter_n) > 1e-10 || Tolerance(K6_iter, K6_iter_n) > 1e-10) {
			K1_iter = K1_iter_n;
			K2_iter = K2_iter_n;
			K3_iter = K3_iter_n;
			K4_iter = K4_iter_n;
			K5_iter = K5_iter_n;
			K6_iter = K6_iter_n;
			
			K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + 2.0 / 3.0 * h, K);
			K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h / 6.0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + 1.0 / 6.0 * h, K);
			K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -15.0 / 8.0 * h)), Scale(K2_iter, 21.0 / 8.0 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + 0.75 * h, K);
			K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -9.0 * h)), Scale(K2_iter, 75.0 / 7.0 * h)), Scale(K3_iter, -5.0 / 7.0 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + h, K);
			K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -3.0 * h)), Scale(K2_iter, 34257.0 / 8750.0 * h)), Scale(K3_iter, -114.0 / 875.0 * h)), Scale(K4_iter, 19.0 / 1250.0 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), time[i - 1] + h * 0.8, K);
			K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 123.0 / 380.0 * h)), Scale(K3_iter, 2.5 * h)), Scale(K4_iter, 3.0 / 20.0 * h)), Scale(K5_iter, -75.0 / 38.0 * h)), Scale(K6_iter, 0)), time[i - 1] + h, K);
		}

		x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Scale(K1_iter_n, 0), Scale(K2_iter_n, 54.0 /133.0)), Scale(K3_iter_n, 32.0 / 21.0)), Scale(K4_iter_n, 1.0 / 18.0)), Scale(K5_iter_n, - 125.0 / 114.0)), Scale(K6_iter_n, 1.0 / 9.0)), h), x_solved[i - 1]);
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}

// R-K Tsitouras 65  Exp..........
vector<vector <double>> ODE_Solve_RKT65(vector<double>(*Fun)(const vector<double>& x0, double t, const vector<double>& K), double t_range[2], double step_size, vector<double> x0, const vector<double>& K) {
	int n = x0.size();
	double h = step_size;
	int len = ceil((t_range[1] - t_range[0]) / h + 1);
	vector<double> time(len);
	vector<vector <double>> x_solved(len, vector<double>(n));
	time[0] = t_range[0];
	x_solved[0] = x0;

	for (int i = 1; i < len; ++i) {
		time[i] = time[i - 1] + h;
		vector<double> K1_iter = Fun(x_solved[i - 1], time[i - 1], K);
		vector<double> K2_iter = Fun(x_solved[i - 1], time[i - 1] + 700.0 / 12757.0 * h, K);
		vector<double> K3_iter = Fun(x_solved[i - 1], time[i - 1] + 487.0 / 5394.0 * h, K);
		vector<double> K4_iter = Fun(x_solved[i - 1], time[i - 1] + 207.0 * h / 2830.0, K);
		vector<double> K5_iter = Fun(x_solved[i - 1], time[i - 1] + h * 4783.0 / 8682.0, K);
		vector<double> K6_iter = Fun(x_solved[i - 1], time[i - 1] + h * 3736.0 / 5943.0, K);
		vector<double> K7_iter = Fun(x_solved[i - 1], time[i - 1] + h * 6631.0 / 6733.0, K);
		vector<double> K8_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);
		vector<double> K9_iter = Fun(x_solved[i - 1], time[i - 1] + h, K);

		vector<double> K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1], K);
		vector<double> K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 700.0 / 12757.0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 700.0 / 12757.0 * h, K);
		vector<double> K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 2277859.0 / 142293924.0 * h)), Scale(K2_iter, 1321930.0 / 17797209.0 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 487.0 / 5394.0 * h, K);
		vector<double> K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 2147197.0 / 125328051.0 * h)), Scale(K2_iter, 23026816.0 / 342342975.0 * h)), Scale(K3_iter, -384743.0 / 34198910.0 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 207.0 * h / 2830.0, K);
		vector<double> K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 29186539.0 / 55002258.0 * h)), Scale(K2_iter, -119621649.0 / 271157072.0 * h)), Scale(K3_iter, 2798806734.0 / 97726297.0 * h)), Scale(K4_iter, -1837543755.0 / 55776281.0 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h * 4783.0 / 8682.0, K);
		vector<double> K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -4016312.0 / 9839773.0 * h)), Scale(K2_iter, 45481347.0 / 269218664.0 * h)), Scale(K3_iter, -1447053163.0 / 71668982.0 * h)), Scale(K4_iter, 3538651296.0 / 143594951.0 * h)), Scale(K5_iter, 34474522.0 / 90990629.0 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h * 3736.0 / 5943.0, K);
		vector<double> K7_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 527404067.0 / 134932851.0)), Scale(K2_iter, 92630627.0 / 108312167.0 * h)), Scale(K3_iter, 1860989051.0 / 82010587.0 * h)), Scale(K4_iter, -1403078476.0 / 52174661.0 * h)), Scale(K5_iter, -138841401.0 / 117689617.0 * h)), Scale(K6_iter, h * 82333737.0 / 51439001.0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h * 6631.0 / 6733.0, K);
		vector<double> K8_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1407304125.0 / 293501414.0 * h)), Scale(K2_iter, 1076954593.0 / 980453270.0 * h)), Scale(K3_iter, 3324345790.0 / 119511403.0 * h)), Scale(K4_iter, -1676391597.0 / 50682665.0 * h)), Scale(K5_iter, -86106570.0 / 53800223.0 * h)), Scale(K6_iter, 203157280.0 / 102205023.0 * h)), Scale(K7_iter, h * 3785594.0 / -184956015.0)), Scale(K8_iter, 0)), time[i - 1] + h, K);
		vector<double> K9_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 18933895.0 / 59901407.0 * h)), Scale(K2_iter, 0)), Scale(K3_iter, 232243646.0 / 77800723.0 * h)), Scale(K4_iter, -3615718163.0 / 1227397738.0 * h)), Scale(K5_iter, 35451175.0 / 164731901.0 * h)), Scale(K6_iter, 73396207.0 / 298398716.0 * h)), Scale(K7_iter, h * 189726995.0 / 268693409.0)), Scale(K8_iter, h * -104915077.0 / 200742721.0)), time[i - 1] + h, K);

		while (Tolerance(K1_iter, K1_iter_n) > 1e-10 || Tolerance(K2_iter, K2_iter_n) > 1e-10 || Tolerance(K3_iter, K3_iter_n) > 1e-10 || Tolerance(K4_iter, K4_iter_n) > 1e-10 || Tolerance(K5_iter, K5_iter_n) > 1e-10 || Tolerance(K6_iter, K6_iter_n) > 1e-10 || Tolerance(K7_iter, K7_iter_n) > 1e-10 || Tolerance(K8_iter, K8_iter_n) > 1e-10 || Tolerance(K9_iter, K9_iter_n) > 1e-10) {
			K1_iter = K1_iter_n;
			K2_iter = K2_iter_n;
			K3_iter = K3_iter_n;
			K4_iter = K4_iter_n;
			K5_iter = K5_iter_n;
			K6_iter = K6_iter_n;
			K7_iter = K7_iter_n;
			K8_iter = K8_iter_n;
			K9_iter = K9_iter_n;

			K1_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1], K);
			K2_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 700.0 / 12757.0)), Scale(K2_iter, 0)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 700.0 / 12757.0 * h, K);
			K3_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 2277859.0 / 142293924.0 * h)), Scale(K2_iter, 1321930.0 / 17797209.0 * h)), Scale(K3_iter, 0)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 487.0 / 5394.0 * h, K);
			K4_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 2147197.0 / 125328051.0 * h)), Scale(K2_iter, 23026816.0 / 342342975.0 * h)), Scale(K3_iter, -384743.0 / 34198910.0 * h)), Scale(K4_iter, 0)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + 207.0 * h / 2830.0, K);
			K5_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 29186539.0 / 55002258.0 * h)), Scale(K2_iter, -119621649.0 / 271157072.0 * h)), Scale(K3_iter, 2798806734.0 / 97726297.0 * h)), Scale(K4_iter, -1837543755.0 / 55776281.0 * h)), Scale(K5_iter, 0)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h * 4783.0 / 8682.0, K);
			K6_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, -4016312.0 / 9839773.0 * h)), Scale(K2_iter, 45481347.0 / 269218664.0 * h)), Scale(K3_iter, -1447053163.0 / 71668982.0 * h)), Scale(K4_iter, 3538651296.0 / 143594951.0 * h)), Scale(K5_iter, 34474522.0 / 90990629.0 * h)), Scale(K6_iter, 0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h * 3736.0 / 5943.0, K);
			K7_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, h * 527404067.0 / 134932851.0)), Scale(K2_iter, 92630627.0 / 108312167.0 * h)), Scale(K3_iter, 1860989051.0 / 82010587.0 * h)), Scale(K4_iter, -1403078476.0 / 52174661.0 * h)), Scale(K5_iter, -138841401.0 / 117689617.0 * h)), Scale(K6_iter, h * 82333737.0 / 51439001.0)), Scale(K7_iter, 0)), Scale(K8_iter, 0)), time[i - 1] + h * 6631.0 / 6733.0, K);
			K8_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 1407304125.0 / 293501414.0 * h)), Scale(K2_iter, 1076954593.0 / 980453270.0 * h)), Scale(K3_iter, 3324345790.0 / 119511403.0 * h)), Scale(K4_iter, -1676391597.0 / 50682665.0 * h)), Scale(K5_iter, -86106570.0 / 53800223.0 * h)), Scale(K6_iter, 203157280.0 / 102205023.0 * h)), Scale(K7_iter, h * 3785594.0 / -184956015.0)), Scale(K8_iter, 0)), time[i - 1] + h, K);
			K9_iter_n = Fun(Add(Add(Add(Add(Add(Add(Add(Add(x_solved[i - 1], Scale(K1_iter, 18933895.0 / 59901407.0 * h)), Scale(K2_iter, 0)), Scale(K3_iter, 232243646.0 / 77800723.0 * h)), Scale(K4_iter, -3615718163.0 / 1227397738.0 * h)), Scale(K5_iter, 35451175.0 / 164731901.0 * h)), Scale(K6_iter, 73396207.0 / 298398716.0 * h)), Scale(K7_iter, h * 189726995.0 / 268693409.0)), Scale(K8_iter, h * -104915077.0 / 200742721.0)), time[i - 1] + h, K);
		}

		x_solved[i] = Add(Scale(Add(Add(Add(Add(Add(Add(Add(Add(Scale(K1_iter_n, -35243815.0 / 177455297.0), Scale(K2_iter_n, 0)), Scale(K3_iter_n, -10187410.0 / 61778981.0)), Scale(K4_iter_n, 3339569.0 / 4973441.0)), Scale(K5_iter_n, 29595575.0 / 112455507.0)), Scale(K6_iter_n, 59599904.0 / 213941517.0 )), Scale(K7_iter_n, 46892093.0 / 154074174.0)), Scale(K8_iter_n, -841.0 / 4121.0)), Scale(K9_iter_n, 1/20.0)), h), x_solved[i - 1]);
	}

	vector<vector <double>> result(n + 1);
	result[0] = time;
	for (int p = 0; p < n; ++p) {
		for (int l = 0; l < len; ++l) {
			result[p + 1].push_back(x_solved[l][p]);
		}
	}

	return result;
}