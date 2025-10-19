// #include <cmath>
// #include "qss/distributions.hpp"
// #include "qss/math_utils.hpp"

// namespace qss{

// // NFL
// double nfl_distribution(double e_k, double mu, double beta) {
// 	double x = beta * (e_k - mu);
// 	if (x > UP_EXP_LIMIT) return 0.0;
// 	if (x < DOWN_EXP_LIMIT) return 1.0;
// 	double expbeta = std::exp(x);
// 	return (expbeta - std::sqrt(expbeta * expbeta + 8.0 * expbeta) + 8.0) / (expbeta + 8.0);
// }

// Vec nfl_distribution(const Vec &e_k, double mu, double beta) {
// 	Vec x; x.reserve(e_k.size());
// 	for (double ek : e_k) x.push_back(beta * (ek - mu));
// 	Vec result = handle_overflow_vector(x, UP_EXP_LIMIT, DOWN_EXP_LIMIT, 1.0, 0.0);
// 	for (std::size_t i = 0; i < x.size(); ++i) {
// 		if (x[i] >= DOWN_EXP_LIMIT && x[i] <= UP_EXP_LIMIT) {
// 			double expbeta = std::exp(x[i]);
// 			result[i] = (expbeta - std::sqrt(expbeta * expbeta + 8.0 * expbeta) + 8.0) / (expbeta + 8.0);
// 		}
// 	}
// 	return result;
// }

// // Implicit NFL_i: we solve the cubic p0*f^3 + p1*f^2 + p2*f + p3 = 0 for a root in [0,1]
// static bool find_root_on_unit_interval(double p0, double p1, double p2, double p3, double &root) {
// 	auto f = [&](double fv) {
// 		return ((p0 * fv + p1) * fv + p2) * fv + p3; // p0*fv^3 + p1*fv^2 + p2*fv + p3
// 	};
// 	double a = 0.0, b = 1.0;
// 	double fa = f(a), fb = f(b);
// 	if (std::isnan(fa) || std::isnan(fb)) return false;
// 	if (fa == 0.0) { root = a; return true; }
// 	if (fb == 0.0) { root = b; return true; }
// 	if (fa * fb > 0.0) return false; // no sign change

// 	// Bisection
// 	for (int iter = 0; iter < 80; ++iter) {
// 		double m = 0.5 * (a + b);
// 		double fm = f(m);
// 		if (std::abs(fm) < 1e-12) { root = m; return true; }
// 		if (fa * fm <= 0.0) { b = m; fb = fm; }
// 		else { a = m; fa = fm; }
// 	}
// 	root = 0.5 * (a + b);
// 	return true;
// }

// double nfl_i_distribution(double e_k, double mu, double beta) {
// 	double x = beta * (e_k - mu);
// 	if (x > UP_EXP_LIMIT) return 0.0;
// 	if (x < DOWN_EXP_LIMIT) return 1.0;
// 	double expbeta = std::exp(x);
// 	// polynomial coefficients for f: p0*f^3 + p1*f^2 + p2*f + p3 = 0
// 	double p0 = expbeta;
// 	double p1 = 2.0 * expbeta + 8.0;
// 	double p2 = -16.0;
// 	double p3 = 8.0;
// 	double froot = 0.0;
// 	if (find_root_on_unit_interval(p0, p1, p2, p3, froot)) {
// 		if (froot >= 0.0 && froot <= 1.0) return froot * 2.0;
// 	}
// 	return 0.0;
// }

// Vec nfl_i_distribution(const Vec &e_k, double mu, double beta) {
// 	Vec x; x.reserve(e_k.size());
// 	for (double ek : e_k) x.push_back(beta * (ek - mu));
// 	Vec result = handle_overflow_vector(x, UP_EXP_LIMIT, DOWN_EXP_LIMIT, 1.0, 0.0);
// 	for (std::size_t i = 0; i < x.size(); ++i) {
// 		if (x[i] >= DOWN_EXP_LIMIT && x[i] <= UP_EXP_LIMIT) {
// 			double expbeta = std::exp(x[i]);
// 			double p0 = expbeta;
// 			double p1 = 2.0 * expbeta + 8.0;
// 			double p2 = -16.0;
// 			double p3 = 8.0;
// 			double froot = 0.0;
// 			if (find_root_on_unit_interval(p0, p1, p2, p3, froot) && froot >= 0.0 && froot <= 1.0)
// 				result[i] = froot * 2.0;
// 			else
// 				result[i] = 0.0;
// 		}
// 	}
// 	return result;
// }
// }