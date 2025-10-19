// #include <cmath>
// #include "qss/distributions.hpp"
// #include "qss/math_utils.hpp"

// namespace qss{

// // FLB
// double flb_distribution(double e_k, double mu, double beta) {
// 	double x = beta * (e_k - mu);
// 	if (x > UP_EXP_LIMIT) return 0.0;
// 	if (x < DOWN_EXP_LIMIT) return 2.0;
// 	return 4.0 / (2.0 + std::exp(x));
// }

// Vec flb_distribution(const Vec &e_k, double mu, double beta) {
// 	Vec x; x.reserve(e_k.size());
// 	for (double ek : e_k) x.push_back(beta * (ek - mu));
// 	Vec result = handle_overflow_vector(x, UP_EXP_LIMIT, DOWN_EXP_LIMIT, 2.0, 0.0);
// 	for (std::size_t i = 0; i < x.size(); ++i) {
// 		if (x[i] >= DOWN_EXP_LIMIT && x[i] <= UP_EXP_LIMIT) result[i] = 4.0 / (2.0 + std::exp(x[i]));
// 	}
// 	return result;
// }

// }

