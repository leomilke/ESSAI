#include "AE_Solve_Algorithm.h"
#include "ODE_Solve_Algorithm.h"
// Newton methods ?
// Homotopy to ODE IVP ?
// Anderson Acceleration ?
// Trust Region ?
// Random, Multi-target in QPSO ?

unsigned int seed = 2024;

// Hybrid Function Homotopy H(x,t)=t*F(x)+(1-t)*G(x) 
vector<double> Fun_Homo(vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), vector<double> x, vector<double> x0, double t, double coef[], vector<double> K) {
	vector<double> Fun_Homo_N = Add(Fun(x, K), Scale(Fun(x0, K), t - 1.0)); // Newton G(x)=F(x)-F(x0)
	vector<double> Fun_Homo_F = Add(Scale(Fun(x, K), t), Scale(Substract(x, x0), 1.0 - t)); // Fixed Point G(x)=x-x0
	vector<double> Fun_Homo_E = Add(Fun(x, K), Scale(Scale(Fun(x0, K), exp(-10.0 * t)), t - 1.0)); // Exponential G(x)=F(x)-exp(-k*t)F(x0)
	return Add(Add(Scale(Fun_Homo_N, coef[0]), Scale(Fun_Homo_F, coef[1])), Scale(Fun_Homo_E, coef[2]));
}

// Function Steffensen
vector<double> Fun_St(vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), vector<double> x, vector<double> x0, double t, double coef[], vector<double> K) {
	vector<double> fai = Add(Fun_Homo(Fun, x, x0, t, coef, K), x); //Function Reform
	return Substract(x, Div(Square(Substract(fai, x)), Add(Substract(Add(Fun_Homo(Fun, fai, x0, t, coef, K), fai), Scale(fai, 2.0)), x)));
}

// Linear Search
vector<double> LS(vector<double> v_start, vector<double> v_end, double num, vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), vector<double> x0, double t, double coef[], vector<double> K) {
	vector<double> step = Scale(Substract(v_end, v_start), 1.0 / num);
	double min = 1e20;
	double index = num;
	for (double i = 1.0; i <= num; ++i) {
		double v = Max(Fun_Homo(Fun, Add(v_start, Scale(step, i)), x0, t, coef, K));
		if (v < min) {
			min = v;
			index = i;
		}
	}
	return Add(v_start, Scale(step, index));
}

// Heuristic algorithm to estimate initial point
// Quasi-PSO (Particle Swarm Optimization)
vector<double> QPSO(double x_abs_max, int n, vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), vector<double> K, double positive, int iter_num, bool fixed_seed) {
	int p_num = 500;
	double alpha = 0.2, d_alp = (0.8 - alpha) / (iter_num - 1.0);
	vector<vector<double>> P(p_num, vector<double>(n));
	
	double min = 1e20, e;
	int index = 0;
	random_device d;
	if (!fixed_seed) { seed = d(); }
	mt19937 G(seed);
	uniform_real_distribution<> Rd(-x_abs_max * (1.0 - double(positive > 0.0)) + positive * double(positive > 0.0), x_abs_max);
	vector<double> P1(n), P2(n), P_diff(n);

	for (int i = 0; i < p_num; ++i) {
		for (int j = 0; j < n; ++j) {
			P[i][j] = Rd(G);
		}
		e = Max(Fun(P[i], K));
		if (min > e) { min = e; index = i; };
	}
	for (; iter_num > 0; --iter_num) {
		for (int i = 0; i < p_num; ++i) {
			if (i != index) {
				P_diff = Substract(P[index], P[i]);
				P1 = Add(P[i], Scale(P_diff, alpha));
				P2 = Add(P[i], Scale(P_diff, 1.0 - alpha));
				if (Max(Fun(P1, K)) > Max(Fun(P2, K))) { P[i] = P2; }
				else { P[i] = P1; }
			}
		}
		for (int i = 0; i < p_num; ++i) {
			e = Max(Fun(P[i], K));
			if (min > e) { min = e; index = i; };
		}
		alpha = alpha + d_alp;
	}

	return P[index];
}

// Adaptive-Step Homotopy Steffensen with Heuristic x0 Estimation & Linear Search
vector<vector<double>> AE_Solve_AHS(vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), double range_abs_max, vector<double> x0, vector<double> K, bool use_PSO) {
	int dim = x0.size();
	
	if (use_PSO) {
		vector<double> x0_est = QPSO(range_abs_max, dim, Fun, K, 0.0, 5, false);
		if (Max(Fun(x0, K)) > Max(Fun(x0_est, K))) { x0 = x0_est; }
	}
	
	double e_use = 1e-10;
	double ls_num = 1;
	double coef[3] = { 0.1,0.6,0.3 };

	double t = 0.0, e = 0.0, dh_init = 0.1;
	vector<double> xtemp = x0;
	vector<vector<double>> re(2);
	double dh_var = 1.0;

	while (t <= 1.0) {
		t = t + dh_var;
		double min_error = 1e20;
		while (Max(Fun_Homo(Fun, x0, x0, t, coef, K)) > e_use) {
			x0 = LS(x0, Fun_St(Fun, x0, x0, t, coef, K), ls_num, Fun, x0, t, coef, K);
			e = Max(Fun_Homo(Fun, x0, x0, t, coef, K));
			if (min_error > e) { min_error = e; }
			else if (e > 1e4 * min_error || min_error == 1.0e20) {
				t = t - dh_var; dh_var = dh_var * 0.5; x0 = xtemp;
				break;
			}
		}
		vector<double> xtemp = x0;
		if (e <= e_use) { dh_var = dh_init; }
	}

	re[0] = x0;
	re[1] = Fun(x0, K);
	return re;
}

// Function Steffensen Basic
vector<double> Fun_St_Basic(vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), vector<double> x, vector<double> K) {
	vector<double> fai = Add(Fun(x, K), x); //Function Reform
	return Substract(x, Div(Square(Substract(fai, x)), Add(Substract(Add(Fun(fai, K), fai), Scale(fai, 2.0)), x)));
}

// Only Steffensen
vector<vector<double>> AE_Solve_Basic(vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), vector<double> x0, vector<double> K) {
	double e_use = 1e-10;
	vector<vector<double>> re(2);
	
	while (Max(Fun(x0, K)) > e_use) {
		x0 = Fun_St_Basic(Fun, x0, K);
	}
	
	re[0] = x0;
	re[1] = Fun(x0, K);
	return re;
}

// Heuristic Quasi-PSO Solver
vector<vector<double>> AE_Solve_PSO(vector<double>(*Fun)(const vector<double>& x, const vector<double>& K), double range_abs_max, int dim, double e, vector<double> K, double pos) {
	vector<double> x;
	int init_iter_num = 15;
	do {
		x = QPSO(range_abs_max, dim, Fun, K, pos, init_iter_num, true);
		++seed;
	} while (Max(Fun(x, K)) > e);
	
	return { x,Fun(x, K) };
}